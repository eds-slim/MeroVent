import sys
print(sys.path)

#sys.path.insert(0, "..")
from pynomo.nomographer import *

from pynomo.nomographer import *
from math import *
from pyx import *

#text.set(text.LatexEngine)
#text.reset()
#text.preamble(r"\RequirePackage{fontspec}")

import pandas
from pygam import LinearGAM, s, f

from pygam.datasets import hepatitis
import matplotlib.pyplot as plt

import locale
locale.setlocale(locale.LC_ALL, '')

multiplier=1
w=9.5
h=18

N_params_GFR={
        'u_min':10.0,
        'u_max':160.0,
        'function':lambda u:-u,
        'title':r'eGFR [ml/min/1.73m2]',
	'title_draw_center':True,
	'title_distance_center': 1.75,
	'title_opposite_tick': False,
        'tick_levels':3,
        'tick_text_levels':2,
        'text_distance_0': .8,
        'tag':'r2'
                }

N_params_GFR_alt={
        'u_min':10.0,
        'u_max':160.0,
        'function':lambda u:-u,
        'tick_levels':0,
        'tick_text_levels':0,
        'tag':'r2'
                }

N_params_LOGPRO={
        'u_min':10**2.4/multiplier,
        'u_max':10**4.4/multiplier,
        'function':lambda u:log10(multiplier*u),
        'title':r'CSF protein [mg/l]',
	'title_draw_center':True,
	'title_distance_center': 1.5,
	'title_opposite_tick': False,
	'title_extra_angle': 180,
	'turn_relative': True,
        'tick_levels':4,
        'tick_text_levels':2,
        'text_distance_0': .8,
        'tick_side':'left',
        'scale_type':'log smart',
	'text_format':r"$%3.0f$ ",
	'text_format_func': lambda s: f"{s:,}",
        'tag':'r1',
        'extra_params': [{
				'u_min': 10001.0/multiplier, 
                              	'u_max': 10**4.4/multiplier,  
                              	'tick_text_levels': 3,  
                              	},
	#			{
	#			'u_min': 1.1*1000/multiplier, 
        #                    	'u_max': 4.0*1000/multiplier,  
        #                    	'tick_text_levels': 3,  
        #                    	}
				] 
                }

N_params_LOGPRO_alt={
        'u_min':10**2.4/multiplier,
        'u_max':10**4.4/multiplier,
        'function':lambda u:log10(multiplier*u),
        'tick_levels':0,
        'tick_text_levels':0,
        'scale_type':'log smart',
        'tag':'r1',
                }

dict = {}
fig, ax = plt.subplots(2,3)
i=0;
for arg in sys.argv[1:len(sys.argv)]:

	dict[arg] = {}
	df=pandas.read_csv('nomogram_data_MIC_' + arg +'.csv')
	dict[arg]={	'arg':arg,
				'df': df,
				'X': df["dose"].values,
				'm': df["m"].values,
				'b': df["b"].values
			}
	dict[arg]['gam_m'] = LinearGAM(s(0, n_splines=25, spline_order=7, constraints='monotonic_dec')).fit(dict[arg]['X'], dict[arg]['m'])

	dict[arg]['pol_m'] = np.poly1d(np.polyfit(dict[arg]['X'], dict[arg]['m'], 10))
	
	dict[arg]['gam_b'] = LinearGAM(s(0, n_splines=25, spline_order=7)).fit(dict[arg]['X'], dict[arg]['b'])

	dict[arg]['pol_b'] = np.poly1d(np.polyfit(dict[arg]['X'], dict[arg]['b'], 10))

	

	ax[i,0].plot(dict[arg]['X'], dict[arg]['m'], label='data')
	ax[i,0].plot(dict[arg]['X'], dict[arg]['gam_m'].predict(dict[arg]['X']), label='monotonic fit')
	ax[i,0].plot(dict[arg]['X'], dict[arg]['pol_m'](dict[arg]['X']), label='pol')
	ax[i,0].legend()

	ax[i,1].plot(dict[arg]['X'], dict[arg]['b'], label='data')
	ax[i,1].plot(dict[arg]['X'], dict[arg]['gam_b'].predict(dict[arg]['X']), label='monotonic fit')
	ax[i,1].plot(dict[arg]['X'], dict[arg]['pol_b'](dict[arg]['X']), label='pol')
	ax[i,1].legend()

	ax[i,2].plot(dict[arg]['m'], dict[arg]['b'], label='data')
	ax[i,2].plot(dict[arg]['gam_m'].predict(dict[arg]['X']), dict[arg]['gam_b'].predict(dict[arg]['X']), label='monotonic fit')
	ax[i,2].legend()

	col=color.cmyk.WildStrawberry if arg == '2' else color.cmyk.Turquoise
	dict[arg]['N_params_3']={
	'tag':'scale'+arg,
        'u_min':1.5,
        'u_max':17.5,
	'base_stop': 17,
        'function_3':lambda u:dict[arg]['pol_m'](u), #gam1.predict(u)[0],
        'function_4':lambda u:-dict[arg]['pol_b'](u), #gam2.predict(u)[0],
        'title':r'Meropenem',
	'title_y_shift': 0.5,
        'tick_levels':3,
        'tick_text_levels':3,
	'tick_side': 'right' if arg=='2' else 'left',
        'scale_type':'linear smart',
        #'title_draw_center':True,
	#'title_distance_center': 1.5,
	#'title_opposite_tick': False,
	'grid_length_0': .5,
        #'grid_length_1': 0.6/4,
        #'grid_length_2': 0.4/4,
        #'grid_length_3': 0.4/4,
        #'grid_length_4': 0.3/4,
	'text_size_0': text.size.scriptsize,
        'text_size_1': text.size.tiny,
        #'text_size_2': text.size.tiny,
        #'text_size_3': text.size.tiny,
        #'text_size_4': text.size.tiny,
        'text_distance_0': .5333,
        #'text_distance_1': 0.8/4,
        #'text_distance_2': 0.6/4,
        #'text_distance_3': 0.6/4,
        #'text_distance_4': 1.0/4,
	'title_distance_center': 1.5,
	'extra_titles': [{'dx':-6,
             			 'dy':-12,
              			'text':r'g/d',
              			'width':5,
              			'pyx_extra_defs':[text.size(-3), trafo.rotate(0)]
              			}],
        'axis_color': col,
	#'title_color': col,
                }
	

	dict[arg]['block_1_params']={
             'block_type':'type_10',
             'width':w,
             'height':h,
             'f1_params':N_params_LOGPRO,
             'f2_params':N_params_GFR,
             'f3_params':dict[arg]['N_params_3'],
	     'isopleth_values':[[2300,70,'y']],
             }

	col=color.cmyk.Turquoise if arg == '2' else color.cmyk.WildStrawberry
	dict[arg]['N_params_3_alt']={
	'tag':'scale'+arg,
        'u_min':1.5*1000/24,
        'u_max':17.5*1000/24,
        'function_3':lambda u:dict[arg]['pol_m'](u/1000*24),  #gam1.predict(u)[0],
        'function_4':lambda u:-dict[arg]['pol_b'](u/1000*24), #gam2.predict(u)[0],
        'title':'',
        'tick_levels':3,
        'tick_text_levels':3,
	'tick_side': 'left' if arg=='2' else 'right',
        'scale_type':'linear smart',
        'title_draw_center':True,
	'title_opposite_tick': True,
	'grid_length_0': .5,
        #'grid_length_1': 0.6/4,
        #'grid_length_2': 0.4/4,
        #'grid_length_3': 0.4/4,
        #'grid_length_4': 0.3/4,
	'text_size_0': text.size.scriptsize,
        'text_size_1': text.size.tiny,
        #'text_size_2': text.size.tiny,
        #'text_size_3': text.size.tiny,
        #'text_size_4': text.size.tiny,
        'text_distance_0': .5333,
        #'text_distance_1': 1.1/4,
        'text_distance_2': 0.15,
        #'text_distance_3': 0.2,
        #'text_distance_4': 1.0/4,
	'title_distance_center': 1.5,
	'extra_titles': [{'dx':-5,
             			 'dy':-12,
              			'text':r'mg/hr',
              			'width':5,
              			'pyx_extra_defs':[text.size(-3), trafo.rotate(0)]
              			}],
        'axis_color': color.cmyk.Turquoise if arg == '2' else color.cmyk.WildStrawberry,
	'title_color': color.cmyk.Turquoise if arg == '2' else color.cmyk.WildStrawberry,
                }
	

	dict[arg]['block_1_params_alt']={
             'block_type':'type_10',
             'width':w,
             'height':h,
             'f1_params':N_params_LOGPRO_alt,
             'f2_params':N_params_GFR_alt,
             'f3_params':dict[arg]['N_params_3_alt'],
	     'isopleth_values':[[2300,70,'y']],
             }

	i=i+1
	

#plt.show()

keys=list(dict.keys())

if '2' in keys:
	dict['2']['N_params_3']['function_3'] = lambda u:dict['2']['pol_m'](u)
	dict['2']['N_params_3']['function_4'] = lambda u:-dict['2']['pol_b'](u)
if '4' in keys:
	dict['4']['N_params_3']['function_3'] = lambda u:dict['4']['pol_m'](u)
	dict['4']['N_params_3']['function_4'] = lambda u:-dict['4']['pol_b'](u)
if '8' in keys:
	dict['8']['N_params_3']['function_3'] = lambda u:dict['8']['pol_m'](u)
	dict['8']['N_params_3']['function_4'] = lambda u:-dict['8']['pol_b'](u)



main_params={
              'filename':'nomograph_MIC_' + ''.join(sys.argv[1:len(sys.argv)]) + '.pdf',
              'paper_height':h,
              'paper_width':w,
              'block_params':[value['block_1_params'] for key, value in dict.items()]+[value['block_1_params_alt'] for key, value in dict.items()],
              'transformations':[('rotate',0.01),('scale paper',)],
              #'title_str':r'$y-m(z)x-b(z)=0$',
		'make_grid': False,
              }

Nomographer(main_params)
