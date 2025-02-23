"""Ce code permet de generer une interface graphique pour le dispositif FPGA de Sitrineo
Il permet d'initialiser (sitrireset), de faire la configuration (sitriconf)
le demarage (sitristart), et l'acquisiion (daqSoC_v3 -RunNumber [RRRR] -NumEventsToRead [EEEE] 
-triggerSW yes -StepTriggerMonitor [MMMM] -JtagInit work -DataSaveLocal yes -Delay 1)
Finalement il permet aussi d'enregistrer et transmettre les resultats sur le PC (scp -r root@192.168.1.12:DataStore/run0001 1) 
afin de pouvoir analyser les donnees avec l'interface de taf.

On mettra un alias Sitrineo qui se connecte en ssh au dispositif et qui lance ce code pour generer l'interface graphique 
(alias sitrineo_gui='python interface_sitrineo.py & ; ssh root@192.168.1.12')


On pourra evntuellement rajouter un bouton TAF qui permet de lancer l'interface graphique de TAF permettant de faire l'analyse 
(taf -run 1 -cfg ./config/sitrineo-m2-B.cfg -guiw)

On rajoute un onglet Help pour donner les directives pour le fonctionnement de la machine (branchement, alimentation, commande manuelle...)
"""

# Update: 10/02/2025 by R.TECHI and D.SENINA

import os
import tkinter as tk
import customtkinter as ctk
import subprocess
from PIL import Image, ImageTk

""" Pour configurer l'interface: """

## Creation fenetre:
fenetre = ctk.CTk() 
fenetre.geometry("500x500")
fenetre.title("Sitrineo")

#Logo sitrineo:
image_path = "logo sitrineo.png"  
image = ctk.CTkImage(light_image=Image.open(image_path), size=(200, 100))

# Ajouter l'image dans un label
image_label = ctk.CTkLabel(fenetre, image=image,text="")
image_label.pack(pady=20)

## Ecriture des commandes dans le terminal:

def test():
    subprocess.run("ls",shell=True) # Fonction test 
    print("test done")

def reset():
    subprocess.run("sitrireset",shell=True)
    print("reset done")

def config():
    subprocess.run("sitriconf",shell=True)
    print("configuartion done")

def start():
    subprocess.run("sitristart",shell=True)
    print("start done")

def popup_window():
    global run_number_value

    # Creer une fenetre pop-up
    popup = ctk.CTkToplevel(fenetre)
    popup.geometry("400x300")
    popup.title("Acquisition settings")

    # Labels et champs de saisie pour les parametres
    label_run_number = ctk.CTkLabel(popup, text="Run Number:")
    label_run_number.pack(pady=5)
    entry_run_number = ctk.CTkEntry(popup)
    entry_run_number.pack(pady=5)

    label_num_events = ctk.CTkLabel(popup, text="Number of Events:")
    label_num_events.pack(pady=5)
    entry_num_events = ctk.CTkEntry(popup)
    entry_num_events.pack(pady=5)

    label_step_trigger = ctk.CTkLabel(popup, text="Step Trigger Monitor:")
    label_step_trigger.pack(pady=5)
    entry_step_trigger = ctk.CTkEntry(popup)
    entry_step_trigger.pack(pady=5)
    
    def close_popup(): # Fermer le pop-up apres execution
        popup.destroy()
    
    def acq(): # Fait l'acquisition
        global run_number_value
        
        run_number = entry_run_number.get()
        num_events = entry_num_events.get()
        step_trigger = entry_step_trigger.get()

        # Sauvegarder la valeur du run_number dans la variable globale
        run_number_value = run_number

        # Construire la commande avec les parametres saisis
        command = "daqSoC_v3 -RunNumber {0} -NumEventsToRead {1} -triggerSW yes -StepTriggerMonitor {2} -JtagInit work -DataSaveLocal yes -Delay 1".format(run_number,num_events,step_trigger)

        # Executer la commande
        subprocess.run(command, shell=True)
        print("Acquisition done")
        close_popup() #Fermeture du popup

    button_execute = ctk.CTkButton(popup, text="Lancer", command=acq) #bouton qui lance l'acquisition
    button_execute.pack(pady=10) 

def save():
    global run_number_value
    #subprocess.run("quit",shell=True)
    commande = "scp -r root@192.168.1.12:DataStore/run00{0} {0}".format(run_number_value)
    subprocess.run(commande, shell=True)
    print(commande) #test

def TAF():
    global run_number_value
    file_path = tk.filedialog.askopenfilename(filetypes=[("Configuration file", "*.cfg")]) 
    file_name = os.path.basename(file_path) #ex: sitrineo-m2-B.cfg (with magnetic field) sitrineo-m2.cfg (without B)
    commande="taf -run {0} -cfg ./config/{1} -guiw".format(run_number_value,file_name)
    #subprocess.run(commande, shell=True)
    print(commande) #test

## Help files: 
def help_readme():
    file=open("./README DAQ.txt","r",encoding="utf-8")
    content = file.read()
    print(content)

def help_code():
    file=open("./README GUI.txt","r",encoding="utf-8")
    content = file.read()
    print(content)

## Buttons:
button_reset = ctk.CTkButton(fenetre, text="reset",command=reset)
button_reset.pack(padx=20, pady=10)
button_config= ctk.CTkButton(fenetre, text="config",command=config)
button_config.pack(padx=20, pady=10)
button_start= ctk.CTkButton(fenetre, text="start",command=start)
button_start.pack(padx=20, pady=10)
button_acq= ctk.CTkButton(fenetre, text="acquisition",command=popup_window)
button_acq.pack(padx=20, pady=10)
button_save= ctk.CTkButton(fenetre, text="Save and export",command=save)
button_save.pack(padx=20, pady=10)
button_TAF= ctk.CTkButton(fenetre, text="TAF",command=TAF)
button_TAF.pack(padx=20, pady=10)

## Help:
menu_bar = tk.Menu(fenetre)
menu_open = tk.Menu(menu_bar, tearoff=0)
menu_open.add_separator()
menu_open.add_command(label="Readme",command=help_readme) # Ouvre le READMI 
menu_open.add_command(label ="More",command=help_code) # Info sur le code
menu_bar.add_cascade(label="HELP",menu=menu_open) 
fenetre.config(menu=menu_bar)

## Main
fenetre.mainloop()
