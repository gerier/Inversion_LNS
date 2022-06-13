# Inversion_LNS


Dans ce dépôt, se trouve l'implémentation de la discrétisation des équations de Navier-Stokes linéarisées (LNS) sous leur forme primitive, c'est-à-dire avec les variable de la densité, de la pression et de la vitesse.

Le code est scindé en plusieurs fichiers. 

* linearised_navier_stokes.py : relatif aux équations LNS
* discretisation.py : librairie pour la discrétisation spatiale et temporelle
* case[...] : cas test avec la version forward et backward de la résolution
