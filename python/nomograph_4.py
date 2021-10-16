import sys
from pynomo.nomographer import *
from math import *
from pyx import *
import pandas
from pygam import LinearGAM, s, f

multiplier=1
w=9.5
h=18

arg=sys.argv[1]
df=pandas.read_csv('nomogram_data_MIC_' + arg +'.csv')
X=df["dose"].values
m=df["m"].values
b=df["b"].values

gam_m = LinearGAM(s(0, n_splines=25, spline_order=7, constraints='monotonic_dec')).fit(X,m)
pol_m = np.poly1d(np.polyfit(X, m, 10))
	
gam_b = LinearGAM(s(0, n_splines=25, spline_order=7)).fit(X,b)
pol_b = np.poly1d(np.polyfit(X,b, 10))


N_params_GFR={
        'tag':'GFR',
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
	}

N_params_GFR_alt={
        'tag':'GFR',
        'u_min':10.0,
        'u_max':160.0,
        'function':lambda u:-u,
        'tick_levels':0,
        'tick_text_levels':0,
}

N_params_LOGPRO={
        'tag':'LOGPRO',
        'u_min':10**2.4/multiplier,
        'u_max':10**4.4/multiplier,
        'function':lambda u:log10(multiplier*u),
        'title':r'CSF protein [mg/l]',
	'title_draw_center':True,
	'title_distance_center': 1.5,
	'title_opposite_tick': False,
	'title_extra_angle': 180,
        'tick_levels':4,
        'tick_text_levels':2,
        'text_distance_0': .8,
        'tick_side':'left',
        'scale_type':'log smart',
	'text_format':r"$%3.0f$ ",
	'text_format_func': lambda s: f"{s:,}",
        'extra_params': [{
		'u_min': 10001.0/multiplier, 
                'u_max': 10**4.4/multiplier,  
                'tick_text_levels': 3,  
                }],
	}

N_params_LOGPRO_alt={
        'tag':'LOGPRO',
        'u_min':10**2.4/multiplier,
        'u_max':10**4.4/multiplier,
        'function':lambda u:log10(multiplier*u),
        'tick_levels':0,
        'tick_text_levels':0,
        'scale_type':'log smart',
	}

col=color.cmyk.WildStrawberry
N_params_3={
	'tag':'scale'+arg,
        'u_min':1.5,
        'u_max':17.5,
	'base_stop': 17,
        'function_3':lambda u:pol_m(u), #gam_m.predict(u)[0],
        'function_4':lambda u:-pol_b(u), #gam_b.predict(u)[0],
        'title':r'Meropenem',
	'title_y_shift': 0.5,
        'tick_levels':3,
        'tick_text_levels':3,
	'tick_side': 'left',
        'scale_type':'linear smart',
        #'title_draw_center':True,
	'title_distance_center': 1.5,
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
	'extra_titles': [{'dx':-6,
		'dy':-12,
              	'text':r'g/d',
              	'width':5,
              	'pyx_extra_defs':[text.size(-3), trafo.rotate(0)]
              	}],
        'axis_color': col,
	}
	
block_1_params={
             'block_type':'type_10',
             'width':w,
             'height':h,
             'f1_params':N_params_LOGPRO,
             'f2_params':N_params_GFR,
             'f3_params':N_params_3,
	     'isopleth_values':[[2300,70,'y']],
             }

col=color.cmyk.MidnightBlue
N_params_3_alt={
	'tag':'scale'+arg,
        'u_min':1.5*1000/24,
        'u_max':17.5*1000/24,
        'function_3':lambda u:pol_m(u/1000*24),  #gam_m.predict(u)[0],
        'function_4':lambda u:-pol_b(u/1000*24), #gam_b.predict(u)[0],
        'title':'',
        'tick_levels':3,
        'tick_text_levels':3,
	'tick_side': 'right',
        'scale_type':'linear smart',
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
	'extra_titles': [{'dx':-5,
		'dy':-12,
              	'text':r'mg/hr',
              	'width':5,
              	'pyx_extra_defs':[text.size(-3), trafo.rotate(0)]
              	}],
        'axis_color': col,
	}
	

block_1_params_alt={
             'block_type':'type_10',
             'width':w,
             'height':h,
             'f1_params':N_params_LOGPRO_alt,
             'f2_params':N_params_GFR_alt,
             'f3_params':N_params_3_alt,
	     'isopleth_values':[[2300,70,'y']],
             }
	

main_params={
	'filename':'nomograph_MIC_' + sys.argv[1] + '.pdf',
        'paper_height':h,
	'paper_width':w,
        'block_params':[block_1_params, block_1_params_alt],
        'transformations':[('rotate',0.01),('scale paper',)],
	}

Nomographer(main_params)
