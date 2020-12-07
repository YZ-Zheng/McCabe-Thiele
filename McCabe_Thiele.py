"""
A Tool to calculate the number of ideal stages using McCabe Thiele graphical method.
"""

import numpy as np
import matplotlib.pyplot as plt
from shapely.geometry import LineString
from scipy.interpolate import make_interp_spline


def diagonal():
    """
    plot the diagonal line y=x

    """
    plt.plot([0, 1], [0, 1], color='black')
    return [0, 1], [0, 1]


def equilibrium_line(x=None, y=None, alpha=0.0):
    """
    Plot the equilibrium line using either x and y values or relative volatility alpha

    """
    # simple input checks
    if x is None:
        x = []
    if y is None:
        y = []
    if x == 0 and alpha == 0:
        print("Invalid input")
        return

    if y == 0 and alpha == 0:
        print("Invalid input")
        return

    # convert to numpy array
    x = np.array(x)
    y = np.array(y)

    # constructing equilibrium line
    if alpha == 0:
        xnew = np.linspace(x.min(), x.max(), 300)
        spl = make_interp_spline(x, y, k=3)
        y_smooth = spl(xnew)
        plt.plot(xnew, y_smooth, color='black')
        return xnew, y_smooth
    else:
        x = np.linspace(0, 1, 300)
        y = np.divide(alpha * x, 1 + (alpha - 1) * x)
        plt.plot(x, y, color='black')
        return x, y


def top_operating_line(xd, rd):
    """
    return the intercept of the TOL

    """
    intercept = xd / (rd + 1)
    return intercept


def number_of_stages(xd, xb, z=0.0, q=0.0, rd=0.0, x=None, y=None, alpha=0.0, total_reflux=False):
    """
    return the number of ideal stages required to achieve the given separation
    """
    if y is None:
        y = []
    if x is None:
        x = []
    plt.plot([xd, xd], [0, xd], color='black')
    plt.plot([xb, xb], [0, xb], color='black')
    plt.margins(x=0, y=0)
    plt.text(xd, 0, "Xd = " + str(xd), fontsize=16)
    plt.text(xb, 0, "Xb = " + str(xb), fontsize=16)
    diagonal()

    if total_reflux and alpha == 0:
        x_equilibrium, y_equilibrium = equilibrium_line(x, y)
        stages = 0
        x = xd

        # calculating No. of stages
        while x > xb:
            line_1 = LineString(np.column_stack((x_equilibrium, y_equilibrium)))
            line_2 = LineString(np.column_stack(([x, 0], [x, x])))
            intersection = line_1.intersection(line_2)
            plt.plot([x, intersection.x], [x, intersection.y], color='orange')
            plt.plot([intersection.x, intersection.x], [intersection.y, intersection.x], color='orange')
            x = intersection.x
            stages += 1

        plt.title("Number of stages = " + str(stages), fontsize=20)

        return stages

    elif total_reflux:
        x_equilibrium, y_equilibrium = equilibrium_line(alpha=alpha)
        stages = 0
        x = xd

        # calculating No. of stages
        while x > xb:
            line_1 = LineString(np.column_stack((x_equilibrium, y_equilibrium)))
            line_2 = LineString(np.column_stack(([x, 0], [x, x])))
            intersection = line_1.intersection(line_2)
            plt.plot([x, intersection.x], [x, intersection.y], color='orange')
            plt.plot([intersection.x, intersection.x], [intersection.y, intersection.x], color='orange')
            x = intersection.x
            stages += 1

        plt.title("Number of stages = " + str(stages), fontsize=20)

        return stages

    elif alpha == 0:
        x_equilibrium, y_equilibrium = equilibrium_line(x, y)

        # constructing the operating lines
        tol_intercept = top_operating_line(xd, rd)
        line_1 = LineString(np.column_stack(([xd, 0], [xd, tol_intercept])))
        line_2 = LineString(np.column_stack(([z, 1], [z, ((1 - z) * q) / (q - 1 + 1e-10) + z])))
        intersection = line_1.intersection(line_2)
        x_tol = [xd, intersection.x]
        y_tol = [xd, intersection.y]
        x_bol = [xb, intersection.x]
        y_bol = [xb, intersection.y]
        plt.plot([xd, intersection.x], [xd, intersection.y], color='blue')
        plt.text(xd, xd, "Top OL", color='blue')
        plt.plot([z, intersection.x], [z, intersection.y], color='blue')
        plt.text(z, z - 0.03, "Feed line", color='blue')
        plt.plot([xb, intersection.x], [xb, intersection.y], color='blue')
        plt.text(xb, xb - 0.026, "Btm OL", color='blue')

        # calculating No. of stages
        stages = 0
        x = xd
        y = xd
        while x > xb:
            if x > x_tol[1]:
                line_1 = LineString(np.column_stack((x_equilibrium, y_equilibrium)))
                line_2 = LineString(np.column_stack(([x, 0], [y, y])))
                intersection = line_1.intersection(line_2)
                x_temp = intersection.x
                y_temp = intersection.y
                plt.plot([x, intersection.x], [y, intersection.y], color='orange')

                line_1 = LineString(np.column_stack((x_tol, y_tol)))
                line_2 = LineString(np.column_stack(([intersection.x, intersection.x], [intersection.y, 0])))
                intersection = line_1.intersection(line_2)

                if not intersection:
                    line_1 = LineString(np.column_stack((x_bol, y_bol)))
                    line_2 = LineString(np.column_stack(([x_temp, x_temp], [y_temp, 0])))
                    intersection = line_1.intersection(line_2)
                    plt.plot([intersection.x, intersection.x], [intersection.y, y_temp], color='orange')
                    x = intersection.x
                    y = intersection.y
                    stages += 1
                else:
                    plt.plot([intersection.x, intersection.x], [intersection.y, y_temp], color='orange')
                    x = intersection.x
                    y = intersection.y
                    stages += 1
            else:
                line_1 = LineString(np.column_stack((x_equilibrium, y_equilibrium)))
                line_2 = LineString(np.column_stack(([x, 0], [y, y])))
                intersection = line_1.intersection(line_2)
                x_temp = intersection.x
                y_temp = intersection.y
                plt.plot([x, intersection.x], [intersection.y, intersection.y], color='orange')

                line_1 = LineString(np.column_stack((x_bol, y_bol)))
                line_2 = LineString(np.column_stack(([intersection.x, intersection.x], [intersection.y, 0])))
                intersection = line_1.intersection(line_2)

                if not intersection:
                    plt.plot([x_temp, x_temp], [x_temp, y_temp], color='orange')
                    stages += 1
                    x = x_temp

                else:
                    plt.plot([intersection.x, intersection.x], [intersection.y, y_temp], color='orange')

                    y = intersection.y
                    x = intersection.x
                    stages += 1

            plt.title("Number of stages = " + str(stages), fontsize=20)

        return stages

    else:
        x_equilibrium, y_equilibrium = equilibrium_line(alpha=alpha)

        # constructing the operating lines
        tol_intercept = top_operating_line(xd, rd)
        line_1 = LineString(np.column_stack(([xd, 0], [xd, tol_intercept])))
        line_2 = LineString(np.column_stack(([z, 1], [z, ((1 - z) * q) / (q - 1 + 1e-10) + z])))
        intersection = line_1.intersection(line_2)
        x_tol = [xd, intersection.x]
        y_tol = [xd, intersection.y]
        x_bol = [xb, intersection.x]
        y_bol = [xb, intersection.y]
        plt.plot([xd, intersection.x], [xd, intersection.y], color='blue')
        plt.text(xd, xd, "Top OL", color='blue')
        plt.plot([z, intersection.x], [z, intersection.y], color='blue')
        plt.text(z, z - 0.03, "Feed line", color='blue')
        plt.plot([xb, intersection.x], [xb, intersection.y], color='blue')
        plt.text(xb, xb - 0.026, "Btm OL", color='blue')

        # calculating No. of stages
        stages = 0
        x = xd
        y = xd
        while x > xb:
            if x > x_tol[1]:
                line_1 = LineString(np.column_stack((x_equilibrium, y_equilibrium)))
                line_2 = LineString(np.column_stack(([x, 0], [y, y])))
                intersection = line_1.intersection(line_2)
                x_temp = intersection.x
                y_temp = intersection.y
                plt.plot([x, intersection.x], [y, intersection.y], color='orange')

                line_1 = LineString(np.column_stack((x_tol, y_tol)))
                line_2 = LineString(np.column_stack(([intersection.x, intersection.x], [intersection.y, 0])))
                intersection = line_1.intersection(line_2)

                if not intersection:
                    line_1 = LineString(np.column_stack((x_bol, y_bol)))
                    line_2 = LineString(np.column_stack(([x_temp, x_temp], [y_temp, 0])))
                    intersection = line_1.intersection(line_2)
                    plt.plot([intersection.x, intersection.x], [intersection.y, y_temp], color='orange')
                    x = intersection.x
                    y = intersection.y
                    stages += 1
                else:
                    plt.plot([intersection.x, intersection.x], [intersection.y, y_temp], color='orange')
                    x = intersection.x
                    y = intersection.y
                    stages += 1
            else:
                line_1 = LineString(np.column_stack((x_equilibrium, y_equilibrium)))
                line_2 = LineString(np.column_stack(([x, 0], [y, y])))
                intersection = line_1.intersection(line_2)
                x_temp = intersection.x
                y_temp = intersection.y
                plt.plot([x, intersection.x], [intersection.y, intersection.y], color='orange')

                line_1 = LineString(np.column_stack((x_bol, y_bol)))
                line_2 = LineString(np.column_stack(([intersection.x, intersection.x], [intersection.y, 0])))
                intersection = line_1.intersection(line_2)

                if not intersection:
                    plt.plot([x_temp, x_temp], [x_temp, y_temp], color='orange')
                    stages += 1
                    x = x_temp

                else:
                    plt.plot([intersection.x, intersection.x], [intersection.y, y_temp], color='orange')

                    y = intersection.y
                    x = intersection.x
                    stages += 1

            plt.title("Number of stages = " + str(stages), fontsize=20)

        return stages


# Input the variables
if __name__ == '__main__':

    distillate_mole_fraction = 0.6816
    bottom_mole_fraction = 0.0493
    feed_mole_fraction = 0.164
    feed_quality = 1.1
    reflux_ratio = 1
    x_equilibrium_data = [0, 0.019, 0.0721, 0.0966, 0.1238, 0.1661, 0.2337, 0.2608, 0.3273, 0.3965, 0.5079, 0.5198, 0.5732, 0.6763, 0.7472, 0.8943]
    y_equilibrium_data = [0, 0.17, 0.3891, 0.4375, 0.4704, 0.5089, 0.5445, 0.558, 0.5826, 0.6122, 0.6564, 0.6599, 0.6841, 0.7385, 0.7815, 0.8943]
    relative_volatility = 0
    total_reflux = False

    number_of_stages(distillate_mole_fraction, bottom_mole_fraction, feed_mole_fraction, q=feed_quality,
                     rd=reflux_ratio, x=x_equilibrium_data, y=y_equilibrium_data, alpha=relative_volatility,
                     total_reflux=total_reflux)

    plt.show()
