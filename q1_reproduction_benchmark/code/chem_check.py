#!/usr/bin/env python3
"""Per-file indel profile from epiBED files (R9.4.1 vs R10.4.1 fingerprint).
Usage: python chem_check.py EPIBED_DIR SAMPLES_TSV OUT_TSV [MAX_READS] [NPROC]"""
import csv, glob, gzip, os, re, sys
from multiprocessing import Pool

TOK = re.compile(r'([A-Za-z])(\d*)')
FULL = re.compile(r'(?:[A-Za-z]\d*)+')
MAX_READS = 200000

def ref_len(s):
    return sum((int(n) if n else 1) for c, n in TOK.findall(s) if c != 'i')

def find_col(lines):
    # 1st choice: column whose decoded reference length equals end - start
    for cols in lines:
        try:
            span = int(cols[2]) - int(cols[1])
        except ValueError:
            continue
        for j in range(3, len(cols)):
            if cols[j] and FULL.fullmatch(cols[j]) and ref_len(cols[j]) == span:
                return j, 'span_match'
    # fallback: longest column that looks like a run-length pattern
    best, blen = None, 0
    for cols in lines:
        for j in range(3, len(cols)):
            if cols[j] and FULL.fullmatch(cols[j]) and len(cols[j]) > blen:
                best, blen = j, len(cols[j])
    return best, 'longest'

def profile(path):
    reads = ref = dbases = devents = ibases = ievents = meth = unmeth = 0
    span_mismatch = 0
    other = {}
    head, pc, how = [], None, ''
    with gzip.open(path, 'rt') as f:
        for line in f:
            if line.startswith('#'):
                continue
            cols = line.rstrip('\n').split('\t')
            if len(cols) < 4:
                continue
            if pc is None:
                head.append(cols)
                if len(head) < 50:
                    continue
                pc, how = find_col(head)
                if pc is None:
                    raise ValueError('no pattern column found; first line: ' + '|'.join(head[0])[:200])
                todo = head
            else:
                todo = [cols]
            for c in todo:
                pat = c[pc] if pc < len(c) else ''
                r = 0
                for code, num in TOK.findall(pat):
                    n = int(num) if num else 1
                    if code == 'i':
                        ibases += n; ievents += 1
                        continue
                    r += n
                    if code == 'd':
                        dbases += n; devents += 1
                    elif code == 'M':
                        meth += n
                    elif code == 'U':
                        unmeth += n
                    elif code not in ('x', 'P'):
                        other[code] = other.get(code, 0) + n
                try:
                    if r != int(c[2]) - int(c[1]):
                        span_mismatch += 1
                except ValueError:
                    pass
                ref += r
                reads += 1
            if reads >= MAX_READS:
                break
    if pc is None and head:          # file shorter than 50 reads
        pc, how = find_col(head)
    kb = ref / 1000.0 if ref else float('nan')
    calls = meth + unmeth
    return {'pattern_col': f'{pc}({how})',
            'reads': reads,
            'mean_span_bp': round(ref / reads, 1) if reads else 'NA',
            'del_bases_per_kb': round(dbases / kb, 2),
            'del_events_per_kb': round(devents / kb, 2),
            'ins_bases_per_kb': round(ibases / kb, 2),
            'ins_events_per_kb': round(ievents / kb, 2),
            'cpg_calls_per_kb': round(calls / kb, 2),
            'frac_meth': round(meth / calls, 4) if calls else 'NA',
            'span_mismatch_reads': span_mismatch,
            'other_codes': ';'.join(f'{k}:{v}' for k, v in sorted(other.items())) or '-'}

def work(path):
    name = os.path.basename(path)[:-len('.epibed.gz')]
    try:
        return name, profile(path)
    except Exception as e:
        return name, {'error': repr(e)}

def main():
    global MAX_READS
    epi_dir, samples_tsv, out_tsv = sys.argv[1:4]
    if len(sys.argv) > 4:
        MAX_READS = int(sys.argv[4])
    nproc = int(sys.argv[5]) if len(sys.argv) > 5 else 16
    meta = {r['sample']: r for r in csv.DictReader(open(samples_tsv), delimiter='\t')}
    files = sorted(glob.glob(os.path.join(epi_dir, '*.epibed.gz')))
    print(f'{len(files)} epiBED files, {MAX_READS} reads each, {nproc} workers', flush=True)
    with Pool(nproc) as pool:
        results = pool.map(work, files)
    fields = ['sample', 'series', 'group', 'batch', 'pattern_col', 'reads', 'mean_span_bp',
              'del_bases_per_kb', 'del_events_per_kb', 'ins_bases_per_kb', 'ins_events_per_kb',
              'cpg_calls_per_kb', 'frac_meth', 'span_mismatch_reads', 'other_codes', 'error']
    rows = []
    for name, res in results:
        m = meta.get(name, {})
        row = {'sample': name, 'series': name.split('.')[0],
               'group': m.get('group', 'other'), 'batch': m.get('batch', '-')}
        row.update(res)
        rows.append(row)
    rows.sort(key=lambda r: (r['group'], r['series'], r['sample']))
    with open(out_tsv, 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=fields, delimiter='\t', restval='')
        w.writeheader()
        w.writerows(rows)
    errs = [r for r in rows if r.get('error')]
    print(f'files with errors: {len(errs)}')
    for r in errs[:5]:
        print('  ', r['sample'], r['error'][:200])
    print('pattern column used:', sorted(set(r.get('pattern_col', '') for r in rows if not r.get('error'))))
    by = {}
    for r in rows:
        if not r.get('error'):
            by.setdefault((r['group'], r['series']), []).append(r)
    print('\nSERIES SUMMARY (medians over files)')
    print('group\tseries\tn_files\tdel_bases/kb\tdel_events/kb\tins_bases/kb\tfrac_meth\tmean_span_bp\tspan_mismatch')
    for (g, s), rs in sorted(by.items()):
        def med(k):
            v = sorted(float(x[k]) for x in rs if x[k] != 'NA')
            return round(v[len(v) // 2], 3) if v else 'NA'
        print(f'{g}\t{s}\t{len(rs)}\t{med("del_bases_per_kb")}\t{med("del_events_per_kb")}\t'
              f'{med("ins_bases_per_kb")}\t{med("frac_meth")}\t{med("mean_span_bp")}\t{med("span_mismatch_reads")}')
    print(f'\nwrote {out_tsv}')

if __name__ == '__main__':
    main()
