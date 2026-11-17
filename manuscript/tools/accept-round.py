#!/usr/bin/env python3
"""Accept a revision round in a marked .tex (E8.4).

Drops the removals (\\textcolor{gray}{...} and {\\color{gray} ...}), unwraps the
insertions (\\textcolor{colR1}{X} -> X, {\\color{colR1} X} -> X), with balanced
braces; comments are kept as is. Nothing else is touched.
\\input files that carry marking get a clean copy with the new suffix, and the
\\input is pointed at it. references_{k} becomes references_{k+1}.

usage: accept-round.py IN.tex OUT.tex [--color colR1] [--from-k 6 --to-k 7]
"""
import argparse, os, re, sys


def match(s, i):
    """Index of the '}' closing the '{' at or after i (comments/escapes skipped)."""
    d = 0; j = i
    while j < len(s):
        c = s[j]
        if c == '\\': j += 2; continue
        if c == '%':
            k = s.find('\n', j); j = len(s) if k < 0 else k; continue
        if c == '{': d += 1
        elif c == '}':
            d -= 1
            if d == 0: return j
        j += 1
    raise ValueError('unbalanced group at offset %d' % i)


def clean(s, color, cnt):
    out = []; i = 0
    tc = '\\textcolor{%s}{' % color
    gc = '{\\color{%s}' % color
    while i < len(s):
        if s[i] == '%' and (i == 0 or s[i-1] != '\\'):
            k = s.find('\n', i); k = len(s) if k < 0 else k
            out.append(s[i:k]); i = k; continue
        if s.startswith('\\color{%s}\n' % color, i) and (i == 0 or s[i-1] == '\n'):
            i += len('\\color{%s}\n' % color); cnt['ins'] += 1; continue  # bare switch inside a float
        if s.startswith('\\textcolor{gray}{', i):
            j = match(s, i + len('\\textcolor{gray}')); i = j + 1; cnt['gray'] += 1; continue
        if s.startswith('{\\color{gray}', i):
            j = match(s, i); i = j + 1; cnt['gray'] += 1; continue
        if s.startswith(tc, i):
            j = match(s, i + len(tc) - 1)
            out.append(clean(s[i+len(tc):j], color, cnt)); i = j + 1; cnt['ins'] += 1; continue
        if s.startswith(gc, i):
            j = match(s, i)
            body = s[i+len(gc):j]
            prev = ''.join(out)[-1:] if out else '\n'
            if body.startswith(' ') and (prev == '' or prev.isspace()):
                body = body[1:]  # the space after \color{..} is ignored at line start
            out.append(clean(body, color, cnt)); i = j + 1; cnt['ins'] += 1; continue
        out.append(s[i]); i += 1
    return ''.join(out)


def marked(s, color):
    return re.search(r'\\textcolor\{(gray|%s)\}|\\color\{(gray|%s)\}' % (color, color),
                     re.sub(r'(?m)(?<!\\)%.*$', '', s)) is not None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('src'); ap.add_argument('dst')
    ap.add_argument('--color', default='colR1')
    ap.add_argument('--from-k', default='6'); ap.add_argument('--to-k', default='7')
    ap.add_argument('--next-color', default='colR2',
                    help='color of the next round, defined after the current one (preamble) and noted in the header')
    ap.add_argument('--date', default='2026-11-16')
    ap.add_argument('--no-inputs', action='store_true',
                    help='only clean the file; do not copy marked \\input files')
    a = ap.parse_args()
    base = os.path.dirname(os.path.abspath(a.src))
    cnt = {'gray': 0, 'ins': 0}
    t = clean(open(a.src, encoding='utf-8').read(), a.color, cnt)
    log = ['%s: gray groups dropped %d, %s groups unwrapped %d'
           % (a.src, cnt['gray'], a.color, cnt['ins'])]
    if not a.no_inputs:
        def sub(m):
            p = m.group(1)
            f = os.path.join(base, p + '.tex')
            if p.startswith('k%s/' % a.to_k) or not os.path.exists(f): return m.group(0)
            s = open(f, encoding='utf-8').read()
            if not marked(s, a.color): return m.group(0)
            q = re.sub(r'_%s$' % a.from_k, '', p) + '_' + a.to_k
            c2 = {'gray': 0, 'ins': 0}
            open(os.path.join(base, q + '.tex'), 'w', encoding='utf-8').write(clean(s, a.color, c2))
            log.append('  %s -> %s: gray %d, %s %d' % (p, q, c2['gray'], a.color, c2['ins']))
            return '\\input{%s}' % q
        t = re.sub(r'(?m)^\\input\{([^}]*)\}', lambda m: sub(m), t)
    if a.next_color and not a.no_inputs:
        d = '\\definecolor{%s}{RGB}{0,90,200}' % a.color
        assert t.count(d) == 1, 'definecolor of %s not found once' % a.color
        t = t.replace(d, d + '\n' + d.replace(a.color, a.next_color))
        t = ('%% k = %s nasce em %s da k = %s com a rodada %s aceita (insercoes mantidas, remocoes apagadas);\n'
             '%% a marcacao seguinte e %s.\n' % (a.to_k, a.date, a.from_k, a.color, a.next_color)) + t
    t = t.replace('\\bibliography{references_%s}' % a.from_k, '\\bibliography{references_%s}' % a.to_k)
    open(a.dst, 'w', encoding='utf-8').write(t)
    print('\n'.join(log))


if __name__ == '__main__':
    main()
