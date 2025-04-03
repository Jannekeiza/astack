
#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt

# Example input data (replace these with actual values)
npoints=100  # Number of points in the trace
sampr=0.01   # Sampling rate
tshft=0.5    # Time shift
swpol=1      # Polarity switch
maxd=1.0e-5  # Maximum amplitude
sname="Trace 1"  # Trace name

# Generate example amplitude data (replace with actual data)
amp = np.random.randn(npoints)

# Apply normalization and DC offset for plotting
if maxd > 1.0e-6:
    amp=swpol * amp / (1.333 * maxd) + np.arange(1, npoints + 1)
else:
    amp=swpol * amp + np.arange(1, npoints + 1)

# Apply time shift
ptim = np.arange(npoints) * sampr + tshft

# Plot trace
plt.figure(figsize=(10, 6))
plt.plot(ptim, amp, label=sname)

# Apply labels to each trace
x = 0.02 * npoints * sampr
z = amp[0] + 0.1  # Adjust label position based on the first amplitude value
plt.text(x, z, sname, fontsize=12, color='red')

# Add plot labels and grid
plt.xlabel("Time (s)")
plt.ylabel("Amplitude")
plt.title("Trace Plot")
plt.grid(True)
plt.legend()

# Show the plot
plt.show()