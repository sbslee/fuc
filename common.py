import copy

class VCFResult():
    def __init__(self):
        self.meta = []
        self.head = ['#CHROM', 'POS', 'ID', 'REF', 'ALT',
                     'QUAL', 'FILTER', 'INFO', 'FORMAT']
        self.data = {}

    @property
    def shape(self):
        return (len(self.head[9:]), len(self.data))

def read_vcf(vcf_path):
    vcf_result = VCFResult()
    with open(vcf_path) as f:
        for line in f:
            if line.startswith('##'):
                vcf_result.meta.append(line)
            elif line.startswith('#CHROM'):
                fields = line.strip().split('\t')
                vcf_result.head = fields
            else:
                fields = line.strip().split('\t')
                record = f'{fields[0]}:{fields[1]}:{fields[3]}:{fields[4]}'
                vcf_result.data[record] = fields
    return vcf_result

def write_vcf(vcf_result, vcf_path):
    with open(vcf_path, 'w') as f:
        f.write(''.join(vcf_result.meta))
        f.write('\t'.join(vcf_result.head) + '\n')
        for record, fields in vcf_result.data.items():
            f.write('\t'.join(fields) + '\n')

def strip_vcf(vcf_result, subfields=None):
    vcf_result2 = VCFResult()
    vcf_result2.head = copy.deepcopy(vcf_result.head)
    for record, fields in vcf_result.data.items():
        fields[2] = '.'
        fields[6] = '.'
        fields[7] = '.'
        if subfields is None:
            _subfields = ['GT']
        else:
            _subfields = ['GT'] + subfields
        index_list = []
        for subfield in _subfields:
            i = fields[8].split(':').index(subfield)
            index_list.append(i)
        fields[8] = ':'.join(_subfields)
        fields[9:] = [':'.join([x.split(':')[i] for i in index_list])
            for x in fields[9:]]
        vcf_result2.data[record] = fields
    return vcf_result2

def merge_vcfs(vcf_result1, vcf_result2, subfields=None):
    vcf1 = strip_vcf(vcf_result1, subfields)
    vcf2 = strip_vcf(vcf_result2, subfields)
    vcf3 = VCFResult()
    vcf3.head = vcf1.head + vcf2.head[9:]

    missing = './.'
    record, fields = next(iter(vcf2.data.items()))
    n = len(fields[8].split(':'))
    for i in range(1, n):
        missing += ':.'

    for record, fields in vcf1.data.items():
        if record in vcf2.data:
            vcf3.data[record] = fields
        else:
            fields = fields + [missing] * vcf2.shape[0]
            vcf3.data[record] = fields
    for record, fields in vcf2.data.items():
        if record in vcf3.data:
            vcf3.data[record] += fields[9:]
        else:
            fields = fields[:9] + [missing] * vcf1.shape[0] + fields[9:]
            vcf3.data[record] = fields
    return vcf3
