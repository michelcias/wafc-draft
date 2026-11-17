#!/usr/bin/env python3
"""Split a clean one-file .tex into a driver plus one file per section (E8.4).

usage: split-sections.py FLAT.tex ms|supp OUTDIR_PARENT [--kdir k7]
  ms   -> ms_7.tex + k7/ms-0-front, ms-1-introduction ... ms-8-back
  supp -> supp_7.tex + k7/supp-0-front, supp-S01-<slug> ...
Cuts happen at lines that start with \\begin{document}, \\section{ (or \\section[) and
\\bibliographystyle; every section file starts at its \\section line.
"""
import os, re, sys

MS_NAMES = ['introduction', 'model', 'theory', 'computation', 'simulation',
            'application', 'discussion']


def slug(title):
    t = re.sub(r'\\texorpdfstring\{([^{}]*(?:\{[^{}]*\}[^{}]*)*)\}\{[^{}]*(?:\{[^{}]*\}[^{}]*)*\}', r'\1', title)
    t = re.sub(r'\\[a-zA-Z]+', ' ', t)
    t = re.sub(r'[^A-Za-z0-9 ]', ' ', t).lower().split()
    stop = {'the', 'of', 'and', 'by', 'a', 'an', 'in'}
    t = [w for w in t if w not in stop][:4]
    return '-'.join(t)


def main():
    flat, kind, outdir = sys.argv[1], sys.argv[2], sys.argv[3]
    kd = 'k7'
    if '--kdir' in sys.argv: kd = sys.argv[sys.argv.index('--kdir') + 1]
    lines = open(flat, encoding='utf-8').read().splitlines(keepends=True)
    ib = next(i for i, l in enumerate(lines) if l.startswith('\\begin{document}'))
    ie = next(i for i, l in enumerate(lines) if l.startswith('\\bibliographystyle'))
    secs = [i for i, l in enumerate(lines) if i > ib and re.match(r'\\section[\[{]', l) and i < ie]
    if kind == 'ms':
        assert len(secs) == 7, len(secs)
        # the back matter starts at the first starred section after the last numbered one
        back = next(i for i, l in enumerate(lines) if i > secs[-1] and l.startswith('\\section*{'))
    else:
        back = ie
    cuts = [ib + 1] + secs + [back]
    names = ['%s-0-front' % kind]
    for n, s in enumerate(secs, 1):
        if kind == 'ms':
            names.append('ms-%d-%s' % (n, MS_NAMES[n - 1]))
        else:
            names.append('supp-S%02d-%s' % (n, slug(re.sub(r'^\\section(\[.*\])?', '', lines[s]))))
    if kind == 'ms':
        names.append('ms-8-back')
    else:
        # supp: the last section runs up to \bibliographystyle, no back file
        cuts = cuts[:-1] + [ie]
    os.makedirs(os.path.join(outdir, kd), exist_ok=True)
    head = ''.join(lines[:ib + 1])
    drv = head
    # ms: front, 7 sections, back (back ends at ie); supp: front, S sections (last ends at ie)
    if kind == 'ms':
        pieces = list(zip(cuts[:-1], cuts[1:])) + [(cuts[-1], ie)]
    else:
        pieces = list(zip(cuts[:-1], cuts[1:]))
    assert len(pieces) == len(names), (len(pieces), len(names))
    for nm, (a, b) in zip(names, pieces):
        open(os.path.join(outdir, kd, nm + '.tex'), 'w', encoding='utf-8').write(''.join(lines[a:b]))
        drv += '\\input{%s/%s}\n' % (kd, nm)
    drv += ''.join(lines[ie:])
    out = 'ms_7.tex' if kind == 'ms' else 'supp_7.tex'
    open(os.path.join(outdir, out), 'w', encoding='utf-8').write(drv)
    print(out, ':', len(names), 'files in', kd)
    for nm in names: print('  ', nm)


if __name__ == '__main__':
    main()
