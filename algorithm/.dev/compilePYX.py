#! /usr/bin/env python
# *-* coding:utf-8 *-*
import os

#os.system('cython _align.pyx')
#raw_input('ok')
os.system('gcc -c -fPIC -I/usr/include/python2.7/ align.c')
os.system('gcc -shared align.o -o align.so')

#raw_input('Successful')
#
#os.system('cython _utils.pyx')
#raw_input('ok')
#os.system('gcc -c -fPIC -I/usr/include/python2.6/ _utils.c')
#os.system('gcc -shared _utils.o -o _utils.so')
#
#raw_input('Successful')
#os.system('cython _cluster.pyx')
#raw_input('ok')
#os.system('gcc -c -fPIC -I/usr/include/python2.6/ _cluster.c')
#os.system('gcc -shared _cluster.o -o _cluster.so')
#
#raw_input('Successful')
