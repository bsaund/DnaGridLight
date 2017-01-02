# import matplotlib.pyplot as plt
import Tkinter as tk
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
    return new_dna

def baseToColor(new_dna):
    dna_color = new_dna
    for base in range(len(new_dna)):
        if new_dna[base] == 'T':
            dna_color[base] = "red"
        elif new_dna[base] == 'C':
            dna_color[base] = "blue"
        elif new_dna[base] == 'A':
            dna_color[base] = "yellow"
        elif new_dna[base] == 'G':
            dna_color[base] = "green"
    return dna_color
    
def weightedChoice(weights,objects):
    return next(compress(objects,multinomial(1, weights,1)[0]))


root = tk.Tk()
canvas = tk.Canvas(root, width=600, height=600, borderwidth=0, highlightthickness=0, bg="black")
canvas.grid()
ind = 0

def _create_circle(self, x, y, r, **kwargs):
    return self.create_oval(x-r, y-r, x+r, y+r, **kwargs)
tk.Canvas.create_circle = _create_circle


def plot_leds(canvas, led_colors):
    """Plots the 2d led_colors array as circles in a grid on a canvas"""
    canvas.delete("all")
    y_off = 100

    for row in led_colors:
        x_off = 100
        for color in row:
            canvas.create_circle(x_off, y_off, 10, fill=color)
            x_off += 50
        y_off += 50

def update():
    """Replots the LEDs circles
    Recalls itself every 100 millis
    """

    global ind
    #leds = [["blue"] * 8 for _ in xrange(10)];
    #leds = [["blue","blue","red","red","blue","blue","blue","blue"] for _ in xrange(10)];
    #leds[0][ind] = "green"
    #ind = (ind + 1) % 8
    new_DNA = HKY85(1,0.3,0.1,.35,.3,.25,.1,DNA)
    color = baseToColor(new_DNA)
    leds = [[color[0],color[1],color[2],color[3]] for _ in xrange(1)];
    plot_leds(canvas, leds)
    root.after(100, update)

#eventually make the time input a variable

mutatedDNA = HKY85(1,0.3,0.1,.35,.3,.25,.1,["A","A","T","A"])
DNA = HKY85(1,0.3,0.1,.35,.3,.25,.1,mutatedDNA)
#color = baseToColor(DNA)
root.after(10, update)    
root.wm_title("LED Grid")
root.mainloop()


