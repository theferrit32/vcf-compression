import sys
import json

def main(args):
    vcf_filename = args[1]
    gt_map = {}
    with open(vcf_filename) as f:
        for line in f:
            line = line.strip()
            if line.startswith('#'):
                continue
            terms = line.split('\t')
            assert(len(terms) > 9)
            for sample_idx in range(9, len(terms)):
                if terms[sample_idx] not in gt_map:
                    gt_map[terms[sample_idx]] = 0
                gt_map[terms[sample_idx]] += 1

    print(json.dumps(gt_map, indent=2))
    with open(vcf_filename + '.gt_freqs.json', 'w') as fout:
        json.dump(gt_map, fout)


if __name__ == '__main__':
    print(sys.argv)
    main(sys.argv)