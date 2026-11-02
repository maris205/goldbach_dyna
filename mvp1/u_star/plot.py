"""Offline AppKit plots; Swift code and compiler cache are confined to u_star."""
from pathlib import Path
import json
import os
import subprocess
import time

HERE=Path(__file__).resolve().parent
OUT=HERE/'results'
WORK=HERE/'plot_work'

SWIFT=r'''
import AppKit
import Foundation
let args=CommandLine.arguments
let s=try JSONSerialization.jsonObject(with:Data(contentsOf:URL(fileURLWithPath:args[1]))) as! [String:Any]
let W=1800.0,H=1100.0
let rep=NSBitmapImageRep(bitmapDataPlanes:nil,pixelsWide:Int(W),pixelsHigh:Int(H),bitsPerSample:8,samplesPerPixel:4,hasAlpha:true,isPlanar:false,colorSpaceName:.deviceRGB,bytesPerRow:0,bitsPerPixel:0)!
NSGraphicsContext.saveGraphicsState();NSGraphicsContext.current=NSGraphicsContext(bitmapImageRep:rep)
NSColor.white.setFill();NSBezierPath(rect:NSRect(x:0,y:0,width:W,height:H)).fill()
func txt(_ t:String,_ x:Double,_ y:Double,_ z:Double=20){(t as NSString).draw(at:NSPoint(x:x,y:y),withAttributes:[.font:NSFont.systemFont(ofSize:z),.foregroundColor:NSColor.black])}
let ymin=s["ymin"] as! Double,ymax=s["ymax"] as! Double
let logy=s["log_y"] as! Bool
func px(_ x:Double)->Double{return 145+(log10(x)-2)/3*1170}
func py(_ y:Double)->Double{let v=logy ? log10(y) : y;return 180+(v-ymin)/(ymax-ymin)*710}
for i in 0...6{let x=145+Double(i)/6*1170;NSColor.lightGray.setStroke();let p=NSBezierPath();p.move(to:NSPoint(x:x,y:180));p.line(to:NSPoint(x:x,y:890));p.stroke();txt(String(format:"%.0f",pow(10,2+Double(i)/2)),x-25,145)}
for i in 0...5{let v=ymin+(ymax-ymin)*Double(i)/5;let y=180+(v-ymin)/(ymax-ymin)*710;NSColor.lightGray.setStroke();let p=NSBezierPath();p.move(to:NSPoint(x:145,y:y));p.line(to:NSPoint(x:1315,y:y));p.stroke();let label=logy ? String(format:"%.0e",pow(10,v)) : String(format:"%.3f",v);txt(label,30,y-8)}
let colors=[NSColor.systemBlue,NSColor.systemOrange,NSColor.systemGreen,NSColor.systemPurple]
let cs=s["curves"] as! [[String:Any]]
for (j,c) in cs.enumerated(){let col=colors[c["color"] as! Int];col.setStroke();col.setFill();let xs=c["x"] as! [Double],ys=c["y"] as! [Double];let p=NSBezierPath();p.lineWidth=2.5
if c["dashed"] as! Bool{p.setLineDash([9,6],count:2,phase:0)}
for i in 0..<xs.count{let point=NSPoint(x:px(xs[i]),y:py(ys[i]));if i==0{p.move(to:point)}else{p.line(to:point)};if !(c["dashed"] as! Bool){NSBezierPath(ovalIn:NSRect(x:point.x-4,y:point.y-4,width:8,height:8)).fill()}}
p.stroke();let ly=860-Double(j)*42;NSBezierPath(rect:NSRect(x:1360,y:ly,width:25,height:4)).fill();txt(c["label"] as! String,1400,ly-8,19)}
if let target=s["target"] as? Double{NSColor.red.setStroke();let p=NSBezierPath();p.lineWidth=2;p.move(to:NSPoint(x:145,y:py(target)));p.line(to:NSPoint(x:1315,y:py(target)));p.stroke();txt("target 0.001",1350,py(target)+4,19)}
txt(s["title"] as! String,145,1025,30);txt(s["ylabel"] as! String,145,945,22);txt("N (logarithmic axis)",560,95,22)
txt(s["note"] as! String,145,45,19);txt(s["note2"] as! String,145,15,19)
NSGraphicsContext.restoreGraphicsState()
try rep.representation(using:.png,properties:[:])!.write(to:URL(fileURLWithPath:args[2]))
'''


def main():
    start=time.perf_counter()
    data=json.loads((OUT/'d4_rerun.json').read_text())
    WORK.mkdir(exist_ok=True)
    swift=WORK/'plot.swift'
    swift.write_text(SWIFT)
    curves=[];errors=[]
    for k in range(1,5):
        rr=[r for r in data['rows'] if r['k']==k]
        for key,label,dashed in [('freq_knead','knead=odo',False),('Haar','Haar',True)]:
            curves.append(dict(label=f'k={k} {label}',color=k-1,dashed=dashed,
                x=[r['N'] for r in rr],y=[r[key] for r in rr]))
        errors.append(dict(label=f'k={k} knead=odo',color=k-1,dashed=False,
            x=[r['N'] for r in rr],y=[r['relative_error_knead'] for r in rr]))
    specs=[('D4_frequency.png',dict(title='Recovered superstable parameter: D4 frequencies',
        ymin=0.,ymax=.55,log_y=False,ylabel='Count / (N - 1)',curves=curves,
        note='k <= 4: supervised finite-period cylinder dictionary; all pairs match addresses.',
        note2='Solid: kneading = odometer. Dashed: Haar. This is not an infinite conjugacy proof.')),
        ('D4_relative_error.png',dict(title='Recovered superstable parameter: relative Haar error',
        ymin=-5.,ymax=0.,log_y=True,ylabel='abs(freq / Haar - 1), logarithmic scale',
        curves=errors,target=.001,
        note='k <= 4: every selected N >= 10000 passes the strict 0.001 target.',
        note2='Maximum at N >= 10000: 0.000900090009 (N=10000, k=3).'))]
    env=dict(os.environ)
    env['CLANG_MODULE_CACHE_PATH']=str(WORK/'cache')
    env['SWIFT_MODULECACHE_PATH']=str(WORK/'cache')
    for name,spec in specs:
        target=OUT/name
        assert not target.exists()
        sp=WORK/(name+'.json')
        sp.write_text(json.dumps(spec))
        subprocess.run(['/usr/bin/swift','-module-cache-path',str(WORK/'cache'),
                        str(swift),str(sp),str(target)],env=env,check=True)
    elapsed=time.perf_counter()-start
    (OUT/'plot_timing.json').write_text(json.dumps(dict(elapsed_seconds=elapsed,backend='Swift/AppKit (no plotting packages downloaded)'),indent=2)+'\n')
    with (HERE/'REPORT.md').open('a') as f:
        f.write('\n## 图表（追加）\n\n[频率与 Haar](results/D4_frequency.png)；[相对误差](results/D4_relative_error.png)。本地无 matplotlib，使用已装 Swift/AppKit 绘制，未下载依赖；编译缓存亦只写本目录。绘图耗时 '+f'{elapsed:.3f}'+' 秒。\n')
    print('Plot seconds',elapsed)


if __name__=='__main__':
    main()
