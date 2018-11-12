# coding: utf-8

import subprocess
process = raw_input ('Please, insert the name of the script you want to execute:\n')
subprocess.Popen(["nohup", "python", str(process)])