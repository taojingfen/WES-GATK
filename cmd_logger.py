#!/usr/bin/env python
# -*- encoding: utf-8 -*-
'''
@File    :   cmd_logger.py
@Time    :   2020/10/10 14:36:01
@Author  :   tao_jingfen
@Contact :   taojingfen@foxmail.com
@Desc    :   Creat a logger file,
             Run the command and output log info into the file
'''

import logging
import subprocess

def create_logger(logger_name: str, logfile: str) -> str:
    """create an instance of logging.

    - Args:
        - logger_name: the name of the logger
        - logfile: the file that write the log info
    - Returns:
        an instance of logging
    """
    logger = logging.getLogger(logger_name)
    logger.setLevel(logging.INFO)
    handler = logging.FileHandler(logfile)
    handler.setLevel(logging.INFO)
    formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
    handler.setFormatter(formatter)
    logger.addHandler(handler)
    return logger

def run_cmd(cmd: str, logger: str, lock=None) -> bool:
    """run the command using subprocess.check_call.

    - Args:
        - cmd: command to run, e.g. ls
        - logger: instance of logging, 
          Any error information will be written to logger
        - lock: multiprocessing.Lock or None
    - Returns:
        True if success, otherwise False
    """
    if lock is not None: lock.acquire()
    logger.info("Command: %s" %cmd)
    if lock is not None: lock.release()

    if lock is not None: lock.acquire()
    cmd_res = subprocess.run(cmd,shell=True,
                        capture_output=True, text=True)
    if lock is not None: lock.release()
    for info in [cmd_res.stderr, cmd_res.stdout]:
        if info != '':
            logger.info(info)
    if cmd_res.returncode == 0:
        return True
    else:
        return False

if __name__ == '__main__':
    logger = create_logger('test','logger_test.txt')
    a = run_cmd('ls -lh',logger)
    print(a)
    logger = create_logger('test2','logger_test.txt')
    a = run_cmd('ls-lh',logger)
    print(a)
