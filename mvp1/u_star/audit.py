"""Read-only SHA-256 audit of the parent tree; writes only below u_star."""
from pathlib import Path
import hashlib
import json
import os
import sys
import time

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
OUT = HERE / 'results'


def digest(p):
    h = hashlib.sha256()
    with p.open('rb') as f:
        for block in iter(lambda: f.read(1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()


def snapshot():
    files, errors = {}, []
    for base, dirs, names in os.walk(ROOT, followlinks=False):
        dirs[:] = sorted(d for d in dirs if Path(base, d) != HERE)
        for name in sorted(names):
            p = Path(base, name)
            if p.is_file():
                try:
                    files[str(p.relative_to(ROOT))] = digest(p)
                except OSError as e:
                    errors.append({'path': str(p), 'error': str(e)})
    return dict(timestamp_unix=time.time(), root=str(ROOT), excluded=str(HERE),
                files=files, errors=errors)


def main():
    OUT.mkdir(exist_ok=True)
    mode = sys.argv[1]
    if mode == 'before':
        data = snapshot()
        target = OUT / 'sha256_before.json'
        assert not target.exists()
    else:
        baseline = json.loads((OUT / 'sha256_before.json').read_text())
        data = snapshot()
        changed = [dict(path=p, before=h, after=data['files'][p])
                   for p, h in baseline['files'].items()
                   if p in data['files'] and data['files'][p] != h]
        logs = [r for r in changed if Path(r['path']).suffix == '.log'
                and 'run' in Path(r['path']).name.lower()]
        data.update(checked=len(baseline['files']),
                    changed_external_run_logs=logs,
                    changed_other=[r for r in changed if r not in logs],
                    missing=sorted(set(baseline['files']) - set(data['files'])),
                    added_outside=sorted(set(data['files']) - set(baseline['files'])),
                    all_existing_unchanged=not changed and not data['errors'])
        target = OUT / ('sha256_' + mode + '.json')
        assert not target.exists()
    with target.open('x') as f:
        json.dump(data, f, ensure_ascii=False, indent=2)
        f.write('\n')
    print(json.dumps({k: v for k, v in data.items() if k != 'files'},
                     ensure_ascii=False, indent=2))
    print('Hashed files:', len(data['files']))


if __name__ == '__main__':
    main()
