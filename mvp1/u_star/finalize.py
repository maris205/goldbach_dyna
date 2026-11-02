"""Append measured completion status; preserve all prior report sections."""
from pathlib import Path
import hashlib
import json
import time

HERE=Path(__file__).resolve().parent
OUT=HERE/'results'


def main():
    root=json.loads((OUT/'u_star.json').read_text())
    cyl=json.loads((OUT/'cylinders.json').read_text())
    d4=json.loads((OUT/'d4_rerun.json').read_text())
    audit=json.loads((OUT/'sha256_after.json').read_text())
    before=json.loads((OUT/'sha256_before.json').read_text())
    timing=json.loads((OUT/'plot_timing.json').read_text())
    verify=json.loads((OUT/'independent_verification.json').read_text())
    success_seconds=root['elapsed_seconds']+d4['elapsed_seconds']+timing['elapsed_seconds']+verify['elapsed_seconds']
    status=dict(u=root['u'],selected_index=root['selected_index'],candidate_count=9,
        identification='Candidate set; main has canonical MS factors RC * RLC * RLLLC * RLLLLLC; no paper-1 identity verification',
        criteria_status=dict(precision_and_residual=True,nine_numerical_candidates=True,
            unique_identification_by_7dp=False,independent_root_and_count_verification=True,
            finite_dictionaries_k1_to_k4=True,dictionaries_k5_to_k7=False,
            D4_primary_large_N_target=True,infinite_conjugacy_proof=False,
            parent_entire_tree_unchanged=audit['all_existing_unchanged']),
        successful_computation_seconds=success_seconds,
        elapsed_wall_seconds_since_baseline=time.time()-before['timestamp_unix'],
        parent_integrity={k:v for k,v in audit.items() if k not in ('files','root','excluded')},
        plots_visually_checked=['D4_frequency.png','D4_relative_error.png'],
        all_task_writes_confined_to=str(HERE),
        no_network_requests=True,no_git_mutations=True)
    report=HERE/'REPORT.md'
    prefix=report.read_bytes()
    with report.open('a') as f:
        f.write('\n## 最终状态与完整性（追加）\n\n')
        f.write('若明确约定 Wₚ=R L^(p−2) C（p≥2），本次主值的 MS 词逐字等于 W₂*W₃*W₅*W₇；这一条件性命名与本地原文缺失是两回事，不能据此断言论文1 的原始 W₇ 约定。只凭七位四舍五入无法唯一定位：六个匹配者均已保留。\n\n')
        f.write('直接以单个符号代替筛标签并不成功。下表在长度 1021020 的周期词上比较两种方向，C 在二值诊断中记为被筛；真正字典采用多字柱词。\n\n')
        f.write('|k|L=存活时匹配率|首错|R=存活时匹配率|首错|\n|---|---|---|---|---|\n')
        for v in cyl['dictionaries']:
            s=v['direct_sign_diagnostic']
            a,b=s['L_as_survive'],s['R_as_survive']
            f.write(f'|{v["k"]}|{a["matching_rate"]:.9f}|{a["first_failure"]}|{b["matching_rate"]:.9f}|{b["first_failure"]}|\n')
        f.write('\nk=5…7 的地址版在当前尺度也不能一概要求 1e-3（尾段相对误差随周期增大）。所选 N≥10000 的最大地址相对误差分别为 '+', '.join(f'{max(r["relative_error_odo"] for r in d4["rows"] if r["k"]==k and r["N"]>=10000):.9g}' for k in (5,6,7))+'；高层拟合符号版相应最大相对误差为 '+', '.join(f'{max(r["relative_error_knead"] for r in d4["rows"] if r["k"]==k and r["N"]>=10000):.9g}' for k in (5,6,7))+'。失败已保留，不扩大 k≤4 的达标结论。\n\n')
        f.write('稳定性量化：七位截断相对主周期词在 n=240 首次出现非临界 L/R 差异；±1e-12 在 n=600 首次出现，±1e-14、±1e-20、±1e-40 在 20001 点内无非临界差异。所有这些非零扰动在 n=210 已与精确临界符号 C 有区别。float64 近似在 100211 点无非临界 L/R 差异，但不能可靠识别 C：套用含 C 的 k=4 主字典仍有 63784 个未知点，首个未知 n=77。因此本次主达标结果使用高精度轨道及明确的临界保护，不把 float64 也称为完整通过。\n\n')
        f.write('初次求根尝试完整算出 9 根后，在“七位四舍五入仅有一个匹配者”的检查处停止；发现实际有六个后改为“匹配者中离转述值最近，且符合给定更细提示”的规则重跑。失败日志 `roots_run.log` 保留；当时未生成结果 JSON，也未向报告追加求根结果。成功求根日志为 `roots_rerun.log`。\n\n')
        f.write('SHA-256 基线共 '+str(audit['checked'])+' 个上级既有文件；缺失 '+str(len(audit['missing']))+'，外部 run 日志变化 '+str(len(audit['changed_external_run_logs']))+'，其他变化 '+str(len(audit['changed_other']))+'。全部既有文件未变='+str(audit['all_existing_unchanged'])+'。变化明细：`'+json.dumps(audit['changed_other'],ensure_ascii=False)+'`；外部 run 日志明细：`'+json.dumps(audit['changed_external_run_logs'],ensure_ascii=False)+'`；运行期间上级新增：`'+json.dumps(audit['added_outside'],ensure_ascii=False)+'`。本次脚本、补丁和命令没有写这些路径；记录为本次工具操作之外的变化，无法由本审计判定具体进程，不覆盖或还原。前后完整清单为 `results/sha256_before.json`、`results/sha256_after.json`；最后可再调用 `audit.py final` 留终态快照。\n\n')
        f.write(f'成功脚本净计算合计 {success_seconds:.3f} 秒；自基线审计至本段追加的墙钟时间 {status["elapsed_wall_seconds_since_baseline"]:.1f} 秒，包含实现、检查与首次求根失败重跑，非纯算法耗时。两张 PNG 已打开检查。所有新代码与产物只写 u_star/；未联网、未 commit/push。\n\n')
        f.write('交付入口：`roots.py`（扫描/求根/星积分解）、`experiment.py`（字典/全量 D4/敏感性）、`verify.py`（Decimal 与 CRT 复核）、`plot.py`、`audit.py`。在上级目录用 `.venv/bin/python -B u_star/<脚本>.py` 运行；结果使用独占创建，已存在时拒绝覆盖，复现实验应先另存本目录再在新的副本中运行。原始数据见 `results/u_star.json`、`results/cylinders.json`、`results/d4_rerun.json`。\n')
    assert report.read_bytes().startswith(prefix)
    with (OUT/'completion.json').open('x') as f:
        json.dump(status,f,ensure_ascii=False,indent=2);f.write('\n')
    outputs=[p for p in HERE.glob('*.py')]+[report]+[p for p in OUT.iterdir() if p.is_file()]
    manifest={str(p.relative_to(HERE)):hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted(outputs)}
    with (OUT/'output_manifest.json').open('x') as f:
        json.dump(manifest,f,ensure_ascii=False,indent=2);f.write('\n')
    print(json.dumps(status,ensure_ascii=False,indent=2))


if __name__=='__main__':
    main()
