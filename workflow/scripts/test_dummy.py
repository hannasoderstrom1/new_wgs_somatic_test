# vim: syntax=python tabstop=4 expandtab
# coding: utf-8

__author__ = "Hanna Soderstrom"
__copyright__ = "Copyright 2023, Hanna Soderstrom"
__email__ = "hanna.soderstrom@gu.se"
__license__ = "GPL-3"


def test_dummy():
    from dummy import dummy
    assert dummy() == 1
