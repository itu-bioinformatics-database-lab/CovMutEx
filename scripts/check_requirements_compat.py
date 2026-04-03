#!/usr/bin/env python3
"""Check pinned requirements against PyPI and suggest versions compatible with a target Python."""
import sys
import json
import urllib.request
from pathlib import Path
try:
    from packaging.version import Version, InvalidVersion
    from packaging.specifiers import SpecifierSet
    _HAS_PACKAGING = True
except Exception:
    # fallback to pkg_resources for version parsing; will not evaluate requires_python
    from pkg_resources import parse_version as _parse_version
    _HAS_PACKAGING = False


def get_pypi_json(package):
    url = f"https://pypi.org/pypi/{package}/json"
    with urllib.request.urlopen(url, timeout=15) as r:
        return json.load(r)


def choose_latest_compatible(releases, target_py):
    """Given releases dict from PyPI, return the latest version string compatible with target_py or None."""
    candidates = []
    for ver, files in releases.items():
        if not files:
            continue
        # If we have packaging available, check requires_python metadata per file when present
        if _HAS_PACKAGING:
            reqs = set()
            for f in files:
                rp = f.get('requires_python')
                if rp:
                    reqs.add(rp)
            if reqs:
                ok = False
                for rp in reqs:
                    try:
                        if SpecifierSet(rp).contains(target_py):
                            ok = True
                            break
                    except Exception:
                        ok = True
                if not ok:
                    continue
            try:
                candidates.append(Version(ver))
            except InvalidVersion:
                continue
        else:
            # packaging not available: accept versions and sort using pkg_resources.parse_version
            try:
                candidates.append(_parse_version(ver))
            except Exception:
                continue
    if not candidates:
        return None
    # sort candidates and return highest
    if _HAS_PACKAGING:
        return str(sorted(candidates)[-1])
    else:
        return str(sorted(candidates, key=lambda v: v)[-1])


def main():
    target_py = sys.argv[1] if len(sys.argv) > 1 else "3.8"
    req_file = Path('genome_extractor/requirements.txt')
    if not req_file.exists():
        print('requirements.txt not found at', req_file)
        sys.exit(2)
    lines = req_file.read_text().splitlines()
    suggestions = []
    for i, line in enumerate(lines, 1):
        s = line.strip()
        if not s or s.startswith('#') or '==' not in s:
            continue
        pkg, pinned = s.split('==', 1)
        pkg = pkg.strip()
        pinned = pinned.strip()
        try:
            data = get_pypi_json(pkg)
        except Exception as e:
            suggestions.append({'line': i, 'package': pkg, 'pinned': pinned, 'error': str(e)})
            continue
        releases = data.get('releases', {})
        # is pinned available at all?
        if pinned in releases and releases.get(pinned):
            # check if pinned has any file with requires_python that disallows target
            files = releases.get(pinned, [])
            disallowed = False
            reqs = set()
            for f in files:
                rp = f.get('requires_python')
                if rp:
                    reqs.add(rp)
            if reqs:
                ok = False
                for rp in reqs:
                    try:
                        if SpecifierSet(rp).contains(target_py):
                            ok = True
                            break
                    except Exception:
                        ok = True
                if not ok:
                    disallowed = True
            if not disallowed:
                continue
        # pinned not present or disallowed — choose alternative
        candidate = choose_latest_compatible(releases, target_py)
        suggestions.append({'line': i, 'package': pkg, 'pinned': pinned, 'suggested': candidate})

    if not suggestions:
        print('All pinned versions appear compatible with Python', target_py)
        return
    print(json.dumps(suggestions, indent=2))


if __name__ == '__main__':
    main()
