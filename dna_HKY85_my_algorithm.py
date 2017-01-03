# import matplotlib.pyplot as plt
import Tkinter as tk
import time
from math import *
import random
from itertools import compress
from numpy.random import multinomial

def HKY85(time,alpha,beta,T,C,A,G,nt):
    h = T + ((T*(A+G))/(T+C))*exp(-beta*time) + (C/(T+C))*exp(-((T+C)*alpha+(A+G)*beta)*time)
    i = C + ((C*(A+G))/(T+C))*exp(-beta*time) - (C/(T+C))*exp(-((T+C)*alpha+(A+G)*beta)*time) 
    j = T + ((T*(A+G))/(T+C))*exp(-beta*time) - (T/(T+C))*exp(-((T+C)*alpha+(A+G)*beta)*time)
    k = C + ((C*(A+G))/(T+C))*exp(-beta*time) + (T/(T+C))*exp(-((T+C)*alpha+(A+G)*beta)*time)
    l = A + ((A*(T+C))/(A+G))*exp(-beta*time) + (G/(A+G))*exp(-((A+G)*alpha+(T+C)*beta)*time)
    m = G + ((G*(T+C))/(A+G))*exp(-beta*time) - (G/(A+G))*exp(-((A+G)*alpha+(T+C)*beta)*time)
    n = A + ((A*(T+C))/(A+G))*exp(-beta*time) - (A/(A+G))*exp(-((A+G)*alpha+(T+C)*beta)*time)
    o = G + ((G*(T+C))/(A+G))*exp(-beta*time) + (A/(A+G))*exp(-((A+G)*alpha+(T+C)*beta)*time)
    p = A * (1-(exp(-beta*time)))
    q = G * (1-(exp(-beta*time)))
    r = T * (1-(exp(-beta*time)))
    s = C * (1-(exp(-beta*time)))

    #probability matrix
    #  T C A G
    #T h i p q
    #C j k p q
    #A r s l m
    #G r s n o

    prob_matrix = [[h,i,p,q],[j,k,p,q],[r,s,l,m],[r,s,n,o]]

    #dna_seq = []
    #for base in range(nt):
    #    base = "ACTG"[random.randint(0,3)]
    #    dna_seq = dna_seq + [base]
    
    #print dna_seq
    dna_seq = nt
    new_dna = dna_seq
    for base in range(len(dna_seq)):
        if dna_seq[base] == 'T':
             new_dna[base] = weightedChoice([h,i,p,q],['T','C','A','G'])
        elif dna_seq[base] == 'C':
             new_dna[base] = weightedChoice([j,k,p,q],['T','C','A','G'])
        elif dna_seq[base] == 'A':
             new_dna[base] = weightedChoice([r,s,l,m],['T','C','A','G'])
        elif dna_seq[base] == 'G':
             new_dna[base] = weightedChoice([r,s,n,o],['T','C','A','G'])
        else:
             print dna_seq
    return new_dna

def baseToColor(dna_seq):
    """Returns the list of colors associated with given list of nucleotides"""
    colormap = {'T':"red", 'C':"blue", 'A':"yellow", 'G':"green"}
    return [colormap[nt] for nt in dna_seq]

def complement(dna_seq):
    """Returns the complement of the given dna sequence"""
    compmap = {'T':'A', 'A':'T', 'C':'G', 'G':'C'}
    return [compmap[nt] for nt in dna_seq]
    
def weightedChoice(weights,objects):
    return next(compress(objects,multinomial(1, weights,1)[0]))


root = tk.Tk()
canvas = tk.Canvas(root, width=600, height=600, borderwidth=0, highlightthickness=0, bg="black")
canvas.grid()

def _create_circle(self, x, y, r, **kwargs):
    return self.create_oval(x-r, y-r, x+r, y+r, **kwargs)
tk.Canvas.create_circle = _create_circle


def plot_leds(canvas, led_colors, comp_colors):
    """Plots the led_colors array as circles in a grid on a canvas"""
    canvas.delete("all")
    y_off = 100
    x_off = 20
    ind = 0;

    for color in led_colors:
        canvas.create_circle(x_off, y_off + 3*sin(ind*.57), 3, fill=color)
        canvas.create_circle(x_off, y_off + 3*sin(ind*.57 + 2.4), 3, fill=comp_colors[ind])
        x_off += 10
        ind += 1

def plot_leds_on_curve(canvas, led_colors, comp_colors, curve):
    """Plots the led_colors array as a dna helix following the given curve
    
    curve -- parametric curve: (x,y) = curve(t)
    tangent -- the derivative of curve: (dx, dy) = tangent(dt)
    """
    canvas.delete("all")
    ind = 0
    eps = 0.0001
    x_off = 300
    y_off = 300
    ratio = 1.0/len(led_colors)

    for ind in xrange(len(led_colors)):
    # for ind in xrange(30,60):
        t = ind*ratio;
        (x,y) = heart_curve(t);
        (dx, dy) = curve(t+eps)

        dx = (dx - x)/eps
        dy = (dy - y)/eps

        dxp = dx/sqrt(dx*dx + dy*dy)
        dy = dy/sqrt(dx*dx + dy*dy)
        dx = dxp
        x = x + x_off
        y = y + y_off
        tx = dy
        ty = -dx


        h = 10*sin(ind*.57 + time.time()*1.5)
        canvas.create_circle(x + tx*h, y + ty*h, 2, fill=led_colors[ind])
        h = 10*sin(ind*.57 + 2.4 + time.time()*1.5)
        canvas.create_circle(x + tx*h, y + ty*h, 2, fill=comp_colors[ind])


def heart_curve(t):
    t = t*2*pi
    x = 16*pow(sin(t),3)
    y = 13*cos(t) - 5*cos(2*t) - 2*cos(3*t) - cos(4*t)
    return (15*x,-15*y)
    
        
def plot_dna(canvas, dna):
    colors = baseToColor(dna)
    comp_colors = baseToColor(complement(dna))
    # plot_leds(canvas, colors, comp_colors)
    plot_leds_on_curve(canvas, colors, comp_colors, heart_curve)

def update():
    """Replots the LEDs circles
    Recalls itself
    """

    global DNA

    DNA = HKY85(2,0.3,0.1,.35,.3,.25,.1,DNA)
    plot_dna(canvas, DNA)
    root.after(30, update)

DNA = HKY85(10,0.3,0.1,.35,.3,.25,.1,["A"]*300)

root.after(10, update)    
root.wm_title("LED Grid")
root.mainloop()


