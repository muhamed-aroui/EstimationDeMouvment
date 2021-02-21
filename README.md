# Estimation De Mouvment 
L'estimation de mouvement consiste à étudier le déplacement des objets dans une séquence
vidéo, en cherchant la corrélation entre deux images successives afin de prédire le changement de position du contenu. Le principe de l’algorithme d’estimation de mouvement 
est de rechercher dans une image référence des similitudes avec l'image source ensuite enregistrer le vecteur de mouvement qui les relie.
Les algorithmes utilisés pour trouver les vecteurs de mouvement se divisent en deux catégories : les méthodes basées sur le pixel et celles basées sur le contenu.

Parmi les méthodes basées sur le pixel, on trouve l’estimation du flux optique caractérisée par
la méthode de Horn et Schunck qui est la première approche de calcul du flot optique, proposé
en 1981 par Horn et Shunck , cette approche consiste à apposer à l‘équation de contrainte du
flux optique un terme global de régularisation pour estimer le champ de vecteurs vitesses. On
procède ensuite à une minimisation à la fois de la contrainte du flux optique et du terme de
régularisation.
## Instruction pour la simulation 

1- ajouter le projet EstimationDeMouvement.cbp dans CodeBlocks

2- simuler le projet

3- vous aurez 4 fichiers .txt 

4- ouvrez le code matlab inclue dans le dossier de projet

5- ajouter les valeurs de U et V ou U_pyr et V_pyr dans les deux matrices U et V dans le code Matlab

6- simuler le code matlab
## Example d'application
### les images d'entrées 
![alt text](https://github.com/muhamed-aroui/EstimationDeMouvment/blob/master/image1.bmp?raw=true)   ![alt text](https://github.com/muhamed-aroui/EstimationDeMouvment/blob/master/image2.bmp?raw=true) 
### Resultats (visualisation via quiver matlab )
![alt text](https://github.com/muhamed-aroui/EstimationDeMouvment/blob/master/MultiResolution.PNG?raw=true)     ![alt text](https://github.com/muhamed-aroui/EstimationDeMouvment/blob/master/zoomedV.PNG?raw=true)
