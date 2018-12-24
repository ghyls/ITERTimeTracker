# ITERTimeTracker
A tool for plotting in 3D or 2D the evolution of a time-dependent IDS, or for visualizing the animation in ParaView. Works with IMAS 3.20.0


# Quick tutorial:

Open "GenerateAll.py". Go to the end of the script.

Edit the options in order to locate the IDS you want to work with. Here is an example


    shot = 22108  
    run = 1000
    user = 'g2mcarpi'                       
    device = 'imas20'           
    version = '3'
    
Run the script with python3, after loading imaenv/3.20.0 (if you are on Gateway)

A "temp_plugin" folder will be generated with the data you need for the plots.

Use the "Track * .py" tools to read the generated files and plot their data (you may need to edit the scripts in otder to select the location of the files).

