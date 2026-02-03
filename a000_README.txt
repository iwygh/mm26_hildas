To reproduce the work in
Matheson and Malhotra 2026
"On the forced orbital plane of the Hilda asteroids":

Run these files in sequential order:
a001_read_nesvorny_elements.py
through
a016_table2.py.
The resulting files will be in 
b001_astorb_horizons.cscv
through
figure5.png.

Because the outputs of integrations from
a004_integrate_real_objects.py 
and
a005_integrate_clones.py 
are very large, I have not included their outputs on Github.
If you wish to reproduce the exact statistics in 
a008_statistics_vmf.py and 
a009_statistics_cf.py,
you can email me to arrange a file exchange.
Otherwise, you can run the integrations yourself and you should get very similar results.

You will need Python 3 with :
Rebound
Numpy
Pandas
Matplotlib
Scipy
Astroquery
Astropy

Iggy Matheson, Monday, February 2, 2026
ianmatheson@arizona.edu