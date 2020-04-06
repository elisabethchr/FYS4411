import operator

import numpy as np
import matplotlib.pyplot as plt

from sklearn.linear_model import LinearRegression
from sklearn.metrics import mean_squared_error, r2_score
from sklearn.preprocessing import PolynomialFeatures

x1 = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0, 2.1, 2.2, 2.3, 2.4]
#y1 = [0.05744014, 0.1581628, 0.05917677, 0.06348687, 0.05255776, 0.06532868, 0.06911291, 0.0202773, 0.0350336, 0.02060009, 0.01558542, 0.01184566, 0.02628105, 0.06427077, 0.01932827, 0.03470512, 0.02133299, 0.01622406, 0.04254373, 0.02178408]
y1 = [0.05744014237627828, 0.15816280351922696, 0.05917676684593566, 0.06348686834839978, 0.052557761433089954, 0.06532867597221427, 0.06911291379105196, 0.020277304720021202, 0.035033595231060366, 0.020600086299230654, 0.015585424673361919, 0.011845657035571992, 0.026281049279027503, 0.06427076843230288, 0.019328265755028867, 0.0347051163795346, 0.021332989955144148, 0.016224057137174015, 0.04254372936547424, 0.021784078105679446, 0.020781258566062236, 0.02314565689589003, 0.02908299340729649, 0.03343456857946587]
x = np.zeros(len(x1))
y = np.zeros(len(x1))

for i in xrange(len(x1)):
    x[i] = x1[i]
    y[i] = y1[i]

print "x = ", x
print "y = ", y

# transforming the data to include another axis
x = x[:, np.newaxis]
y = y[:, np.newaxis]

print "x = ", x
print "y = ", y

polynomial_features= PolynomialFeatures(degree=2)
x_poly = polynomial_features.fit_transform(x)

model = LinearRegression()
model.fit(x_poly, y)
y_poly_pred = model.predict(x_poly)

rmse = np.sqrt(mean_squared_error(y,y_poly_pred))
r2 = r2_score(y,y_poly_pred)
print(rmse)
print(r2)

# plt.plot(x, y, 'b.')
# # sort the values of x before line plot
# sort_axis = operator.itemgetter(0)
# sorted_zip = sorted(zip(x,y_poly_pred), key=sort_axis)
# x, y_poly_pred = zip(*sorted_zip)
# plt.plot(x, y_poly_pred, color='m')
# plt.show()

plt.plot(x, y, 'b.')
# sort the values of x before line plot
sort_axis = operator.itemgetter(0)
sorted_zip = sorted(zip(x,y_poly_pred), key=sort_axis)
x, y_poly_pred = zip(*sorted_zip)
plt.plot(x, y_poly_pred, 'r--')
plt.xlabel(r'$\Delta x$')
plt.ylabel(r'$\sigma_b$')
plt.title('Blocking error vs. step length')
plt.grid('on')
plt.tight_layout()
plt.savefig('../data/plots/d/blocking_error_vs_steplength_fitpol2.pdf')
plt.show()
