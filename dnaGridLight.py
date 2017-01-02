# import matplotlib.pyplot as plt
import Tkinter as tk



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
    leds = [["blue"] * 8 for _ in xrange(10)];
    leds[0][ind] = "green"
    ind = (ind + 1) % 8

    plot_leds(canvas, leds)
    root.after(100, update)


root.after(10, update)    
root.wm_title("LED Grid")
root.mainloop()


