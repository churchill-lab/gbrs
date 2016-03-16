import emase


def hybridize(**kwargs):
    pass


def align(**kwargs):
    pass


def bam2emase(**kwargs):
    pass


def mask(**kwargs):
    pass


def quantify(**kwargs):
    apply_genotype = kwargs.get('apply_genotype')
    if apply_genotype:
        alnmat = mask(alnmat, **kwargs)
    pass
