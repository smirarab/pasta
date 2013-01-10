#! /usr/bin/env python

class TaxaLabelsMismatchError(Exception):
    def __init__(self, *args, **kwargs):
        Exception.__init__(self, *args, **kwargs)

