# Parameters
i = 0
s = 5

# Configure Axis
set size ratio 1
set xrange[0:12.5173]
set yrange[0:12.5173]
set xlabel "X"
set ylabel "Y"

# Setup Color palette
set cbrange [1.2:-0.]
set pal viridis
set tics nomirror

# Plot
plot '../output/velocity.txt' index i u 1:2:3 with image, '../output/velocity.txt' index i u 1:2:($4/(s*sqrt(($4)**2+($5)**2))):($5/(s*sqrt(($4)**2+($5)**2))) every 8:8 with vectors lc -1 filled notitle
