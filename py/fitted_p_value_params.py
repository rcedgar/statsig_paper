import math

name_ref2x1 = {}
name_ref2m0 = {}
name_ref2c0 = {}
name_ref2m = {}
name_ref2c = {}

def set_pvalue_params(name, ref, x1, m0, c0, m, c):
	name_ref2x1[(name, ref)] = x1
	name_ref2m0[(name, ref)] = m0
	name_ref2c0[(name, ref)] = c0
	name_ref2m[(name, ref)] = m
	name_ref2c[(name, ref)] = c

def get_pvalue(name, ref, score):
	x1 = name_ref2x1[(name, ref)]
	if score < x1:
		m = name_ref2m0[(name, ref)]
		c = name_ref2c0[(name, ref)]
	else:
		m = name_ref2m[(name, ref)]
		c = name_ref2c[(name, ref)]
	p = math.exp(m*score + c)
	assert p >= 0
	if p > 1:
		p = 1
	return p

def get_evalue(name, score, dbsize):
	return get_pvalue(name, score)*dbsize

set_pvalue_params(name='tm', ref='fold', x1=0.55, m0=-28, c0=7.1, m=-35, c=11)
set_pvalue_params(name='tm', ref='sf', x1=0.6, m0=-22, c0=5.4, m=-33, c=12)
set_pvalue_params(name='tm', ref='sfc', x1=0.6, m0=-23, c0=5.7, m=-42, c=17)
set_pvalue_params(name='bits', ref='fold', x1=80, m0=-0.097, c0=-2.7, m=-0.021, c=-8.7)
set_pvalue_params(name='bits', ref='sf', x1=1e+02, m0=-0.06, c0=-3.4, m=-0.022, c=-7.2)
set_pvalue_params(name='bits', ref='sfc', x1=1e+02, m0=-0.073, c0=-3.2, m=-0.049, c=-5.6)
set_pvalue_params(name='ts', ref='fold', x1=0.15, m0=-93, c0=1.9, m=-32, c=-7.2)
set_pvalue_params(name='ts', ref='sf', x1=0.15, m0=-78, c0=0.83, m=-30, c=-6.3)
set_pvalue_params(name='ts', ref='sfc', x1=0.11, m0=-80, c0=-0.58, m=-52, c=-3.7)
set_pvalue_params(name='dali', ref='fold', x1=10, m0=-0.86, c0=-1.8, m=-0.38, c=-6.6)
set_pvalue_params(name='dali', ref='sf', x1=10, m0=-0.58, c0=-2.4, m=-0.38, c=-4.3)
set_pvalue_params(name='dali', ref='sfc', x1=10, m0=-0.63, c0=-2.3, m=-0.6, c=-2.5)
