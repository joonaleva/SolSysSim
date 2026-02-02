import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
#import os ; os.system("clear")
from matplotlib.widgets import Slider,Button,TextBox

data = np.loadtxt("output.dat")
#data = np.loadtxt("Simulation_output_file_Pluto_real.dat")

t = np.array(data[:,0])
x = np.array(data[:,1::3])
y = np.array(data[:,2::3])
z = np.array(data[:,3::3])

#distance = []
#for i in range(len(t)):
#    distance.append(np.sqrt(x[i,0]**2+y[i,0]**2+z[i,0]**2))
#print(f'Mean dist = {sum(distance)/len(distance)}')
#print(f'Distance between 1st and last = {np.sqrt((x[0,1]-x[-1,1])**2+(y[0,1]-y[-1,1])**2+(z[0,1]-z[-1,1])**2)}')

ColorIndex = False
color = [                                       # LIST OF PREDEFINED MATPLOTLIB COLORS FOR SOME SOLAR SYSTEM TARGETS
    # name, colour
    ("Sun", "orange"),
    ("Mercury", "slategrey"),
    ("Venus", "navajowhite"),
    ("Earth", "green"),
    ("Mars", "darkorange"),
    ("Jupiter", "sandybrown"),
    ("Saturn", "wheat"),
    ("Uranus", "powderblue"),
    ("Neptune", "royalblue"),
    ("Ceres", "silver"),
    ("Pluto", "tan"),
    ("Eris","slategray"),
    ("Sedna", "midnightblue"),
    ("Oumuamua", "black"),
    ("Moon", "darkgray"),
    ("Io", "khaki"),
    ("Europa", "burlywood"),
    ("Ganymede", "darkgray"),
    ("Callisto", "dimgray"),
    ("Enceladus", "cyan"),
    ("Titan", "gold"),
    ("Triton", "gainsboro"),
    ("Charon", "peru"),
    ("Dysnomia", "darkslategray")
]

targets = []                                    # Read target sequence
try:                                            # Failsafe if plot used without target list
    with open("sequence.txt") as seq:
        for obj in seq:
            obj = obj.rstrip('\n')
            for cname, ccolor in color:
                if obj == cname:
                    ColorIndex = True
                    break
                else: ColorIndex = False
            if ColorIndex:
                targets.append((obj, ccolor))
            else:
                targets.append((obj, "darkgray"))
except FileNotFoundError:
    for i in range(np.size(x, 1)):
        targets.append(("", "darkgray"))

try:                                              # Define units
    with open("config.txt", mode='r') as config:
        txt = config.readline()
        unitInd = txt.split("\t")[1]
        if unitInd == "1": units = ("m", "s", "kg")
        elif unitInd == "3": units = ("pc", "year", "M_sol")
        else: units = ("au", "days", "kg")
except FileNotFoundError:
    units = ("au", "days", "kg")

def SolSys2D():     # 2D plot of the input data
    plt.figure()
    plt.xlabel(f'x ({units[0]})')
    plt.ylabel(f'y ({units[0]})')
    plt.xlim(-50,50)
    plt.ylim(-50,50)

    for i in range(0,len(targets)):
        plt.plot(x[:,i],y[:,i], label=targets[i][0], color=targets[i][1])


    
    plt.legend(shadow=True,fontsize=11)

def SolSys3D():     # 3D plot of the input data
    fig3D = plt.figure()
    ax3D = fig3D.add_subplot(projection='3d')
    ax3D.set_xlabel(f'x ({units[0]})')
    ax3D.set_ylabel(f'y ({units[0]})')
    ax3D.set_zlabel(f'z ({units[0]})')
    ax3D.set_xlim(-50,50)   #50 is about the whole solar systen
    ax3D.set_ylim(-50,50)
    ax3D.set_zlim(-50,50)

    for i in range(0,len(targets)):
        plt.plot(x[:,i],y[:,i], z[:,i], label=targets[i][0], color=targets[i][1])
    
    ax3D.legend(shadow=True,fontsize=11)

def SolPos():       # suns position coordinates as a function of time
    plt.figure()
    plt.title('Solar position')
    plt.xlabel(f't ({units[1]})')
    plt.ylabel(f'pos ({units[0]})')

    plt.plot(t, x[:,0], label="x")
    plt.plot(t, y[:,0], label="y")
    plt.plot(t, z[:,0], label="z")
    plt.legend(shadow=True,fontsize=11)

### ANIMATION BEGINS ###

fig = plt.figure()
ax = fig.add_subplot(projection='3d')

xdata = [None]*len(targets)
ydata = [None]*len(targets)
zdata = [None]*len(targets)

ln=[None]*len(targets)  # Object trails
pos=[None]*len(targets) # Object points

for i in range(0,len(targets)):
    ln[i], = ax.plot([],[],[], '-', color=targets[i][1], alpha=0.5)
    pos[i], = ax.plot([],[],[],'o', label=targets[i][0], color=targets[i][1])
    xdata[i], ydata[i], zdata[i] = [],[],[]

change = False      #variable to check if something is changed so doesnt have to update without changes

# SPEED SLIDER
speed = len(t)/t[-1]
ax_speed = fig.add_axes([0.3, 0.05, 0.45, 0.03])
speed_slider = Slider(ax=ax_speed, label='Speed', valmin=-3, valmax=3, orientation='horizontal', valinit=0)
def on_speed_change(val):
    global speed
    speed = len(t)/t[-1] * 10**val
speed_slider.on_changed(on_speed_change)

# SCALE SLIDER
scale = 50
ax_scale = fig.add_axes([0.1, 0.2, 0.03, 0.6])
scale_slider = Slider(ax=ax_scale, label='Scale', valmin=-4, valmax=4, orientation='vertical', valinit=1.69897)
def on_scale_change(val):
    global scale, change
    scale = 10**val
    change = True
scale_slider.on_changed(on_scale_change)

# PAUSE BUTTON
pause = False
ax_pause = fig.add_axes([0.81, 0.05, 0.1, 0.05])
pause_button = Button(ax_pause, 'Pause')
def on_pause(val):
    global pause, change
    pause = not pause
    change = True
    if pause: pause_button.label.set_text("Resume")
    else: pause_button.label.set_text("Pause") 
pause_button.on_clicked(on_pause)

# RESET BUTTON
ax_reset = fig.add_axes([0.1, 0.05, 0.1, 0.05])
reset_button = Button(ax_reset, 'Reset')
def on_reset(val):
    global scale, frame_index
    #scale = 50
    frame_index = 0
    for i in range(0,len(targets)):
        xdata[i], ydata[i], zdata[i] = [],[],[]
reset_button.on_clicked(on_reset)

# RESET TRAILS BUTTON
ax_reset_trail = fig.add_axes([0.81, 0.15, 0.1, 0.05])
reset_trail_button = Button(ax_reset_trail, 'Reset Trails')
def reset_trail(val):
    for i in range(0,len(targets)):
        xdata[i], ydata[i], zdata[i] = [],[],[]
reset_trail_button.on_clicked(reset_trail)

# HIDE TRAILS BUTTON
Alpha = False
ax_hide_trail = fig.add_axes([0.81, 0.1, 0.1, 0.05])
hide_trail_button = Button(ax_hide_trail, 'Hide Trails')
def on_hide_trail(val):
    global Alpha, change
    Alpha = not Alpha
    change = True
    if Alpha: hide_trail_button.label.set_text("Show Trails")
    else: hide_trail_button.label.set_text("Hide Trails") 
hide_trail_button.on_clicked(on_hide_trail)

# TOP DOWN VIEW BUTTON
ax_dim = fig.add_axes([0.81, 0.9, 0.1, 0.05])
dim_button = Button(ax_dim, 'Top Down View')
def toggle_dim(event):
    global change
    change = True
    ax.view_init(elev=90, azim=-90)     # top-down view
dim_button.on_clicked(toggle_dim)

# TOGGLE LEGEND BUTTON
legend = True
ax_leg = fig.add_axes([0.81, 0.8, 0.1, 0.05])
leg_button = Button(ax_leg, 'Toggle Legend')
def toggle_leg(event):
    global legend, change
    if legend:
        ax.get_legend().remove()
    else:
        ax.legend(shadow=True,fontsize=11)
    legend = not legend
    change = True
leg_button.on_clicked(toggle_leg)

# PLANET TEXT BOX
planet_index = 0
ax_pl_text = fig.add_axes([0.472, 0.94, 0.1, 0.05])
text_box = TextBox(ax_pl_text, label="", textalignment="center")
text_box.set_val("none")
text_box.stop_typing()

# NEXT PLANET VIEW BUTTON
ax_next = fig.add_axes([0.572, 0.94, 0.1, 0.05])
next_button = Button(ax_next, 'Next')
def next_planet(event):
    global planet_index, change
    if planet_index == len(targets): 
        planet_index = 0
        change = True
        text_box.set_val("none")
    else: 
        planet_index += 1
        text_box.set_val(targets[planet_index-1][0])
    text_box.stop_typing()
    if center: reset_trail(0) 
next_button.on_clicked(next_planet)

# PREVIOUS PLANET VIEW BUTTON
ax_prev = fig.add_axes([0.372, 0.94, 0.1, 0.05])
prev_button = Button(ax_prev, 'Previous')
def prev_planet(event):
    global planet_index, change
    if planet_index == 0: 
        planet_index = len(targets)
        text_box.set_val(targets[planet_index-1][0])
    else: 
        planet_index -= 1
        if planet_index==0: 
            change = True
            text_box.set_val("none")
        else: 
            text_box.set_val(targets[planet_index-1][0])
    text_box.stop_typing()
    if center: reset_trail(0)
prev_button.on_clicked(prev_planet)

# CENTER ON PLANET BUTTON
center = False
ax_center = fig.add_axes([0.2, 0.94, 0.1, 0.05])
center_button = Button(ax_center, 'Center on Target')
def center_planet(event):
    global center, change
    center = not center
    change = True
    if center: center_button.label.set_text("Follow Target")
    else: center_button.label.set_text("Center on Target") 
    reset_trail(0)
center_button.on_clicked(center_planet)


# ANIMATION INITIALIZATION
def init():
    global frame_index, scale
    frame_index = 0
    ax.set_xlim(-scale, scale)
    ax.set_ylim(-scale, scale)
    ax.set_zlim(-scale, scale)
    ax.set_xlabel(f"x ({units[0]})")
    ax.set_ylabel(f"y ({units[0]})")
    ax.set_zlabel(f"z ({units[0]})")
    ax.legend(shadow=True,fontsize=11)
    return ln, 

# ANIMATION UPDATE FUNCTION
def update(frame):
    global frame_index, change
    if not pause: frame_index += speed
    if int(frame_index) > len(t): 
        frame_index -= len(t)
        reset_trail(0)
    frame = int(frame_index) % len(t)
    if planet_index == 0:
        for i in range(0,len(targets)):
            xdata[i].append(x[frame,i])
            ydata[i].append(y[frame,i])
            zdata[i].append(z[frame,i])
            ln[i].set_data(xdata[i], ydata[i])
            ln[i].set_3d_properties(zdata[i])
            pos[i].set_data([x[frame,i]], [y[frame,i]])
            pos[i].set_3d_properties([z[frame,i]])
    elif center:
        for i in range(0,len(targets)):
            xdata[i].append(x[frame,i]-x[frame,planet_index-1])
            ydata[i].append(y[frame,i]-y[frame,planet_index-1])
            zdata[i].append(z[frame,i]-z[frame,planet_index-1])
            ln[i].set_data(xdata[i], ydata[i])
            ln[i].set_3d_properties(zdata[i])
            pos[i].set_data([xdata[i][-1]], [ydata[i][-1]])
            pos[i].set_3d_properties([zdata[i][-1]])
    else:
        ax.set_xlim(-scale+x[frame,planet_index-1], scale+x[frame,planet_index-1])
        ax.set_ylim(-scale+y[frame,planet_index-1], scale+y[frame,planet_index-1])
        ax.set_zlim(-scale+z[frame,planet_index-1], scale+z[frame,planet_index-1])

        for i in range(0,len(targets)):
            xdata[i].append(x[frame,i])
            ydata[i].append(y[frame,i])
            zdata[i].append(z[frame,i])
            ln[i].set_data(xdata[i], ydata[i])
            ln[i].set_3d_properties(zdata[i])
            pos[i].set_data([x[frame,i]], [y[frame,i]])
            pos[i].set_3d_properties([z[frame,i]])


    ax.set_title(f'Time: {t[frame]:.3f} {units[1]}')

    if change:
        for i in range(0,len(targets)): 
            if Alpha: 
                ln[i].set_alpha(0)
            else: 
                ln[i].set_alpha(0.5)
        if center or planet_index==0:
            ax.set_xlim(-scale, scale)
            ax.set_ylim(-scale, scale)
            ax.set_zlim(-scale, scale)
        else:
            ax.set_xlim(-scale+x[frame,planet_index-1], scale+x[frame,planet_index-1])
            ax.set_ylim(-scale+y[frame,planet_index-1], scale+y[frame,planet_index-1])
            ax.set_zlim(-scale+z[frame,planet_index-1], scale+z[frame,planet_index-1])

        change = False
    return ln, 

ani = FuncAnimation(fig, update, frames=10000000, init_func=init, blit=False, repeat=True, interval=0)

#SolSys2D()
#SolSys3D()
#SolPos()
plt.show()