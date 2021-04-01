* Se connecter aux serveurs du CSTB sur une machine "de calcul" (tout sauf ena & biplan)
* Créer un répertoire de travail
* Créer un conda env avec le fichier de config et la commande suivante : `conda env create -n hail_env -f hail_env.yml`
* Activer l'environnement : `conda activate hail_env`
* Charger le module comportant java 1.8 (pour spark) : `module load module load java/java-1.8`
* Lancer jupyter lab à partir de l'environnement conda créé (python 3.6 comportant hail) : `jupyter lab` 
* Depuis sa machine, ouvrir un tunnel (Sarah pourra t'expliquer :) ) : 
    * Windows : modifier fichier config SSH : %USERPROFILE%\.ssh\.config 
    * Ajouter la config suivante (exemple avec bipbip, à remplacer par une autre machine si nécessaire): 
        ```
        Host bipbip
          ForwardX11 yes
          HostName bipbip
          User user
          ProxyCommand C:\Windows\System32\OpenSSH\ssh.exe -w 120 -q -W %h:%p user@ssh.lbgi.fr
          IdentitiesOnly true
          LocalForward 8888 localhost:8888
        ```
        
    * Ouvrir un navigateur en local et taper `localhost:8888`
    * Normalement c'est bon, on lancer le notebook `hail_exomes.ipynb`
