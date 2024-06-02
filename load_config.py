#!/usr/bin/env python
# -*- encoding: utf-8 -*-
'''
@File    :   load_config.py
@Time    :   2020/10/10 14:36:01
@Author  :   tao_jingfen
@Contact :   zb-taojingfen@kingmed.com.cn
@Desc    :   get the config information from config.ini
'''

import configparser

def get_config(configfile: str, section: str) -> dict:
    """get the config information of different section in the config file.

    - Args:
        - configfile: the filename of the configfile, e.g. config.ini:
        [test1]
        a=1
        b=2
        [test2]
        c=3
        - section: the section in the file and try to get, e.g. test1
    - Returns:
        the dict of the getting section
        e.g. {'a':'1','b':'2'}
    """
    config = configparser.ConfigParser()
    config.read(configfile)
    return dict(config.items(section))

if __name__ == '__main__':
    database = get_config('wes_config.ini','database')
    print(database)