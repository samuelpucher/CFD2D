# Parameters
start = 0
end = 100
s = 5

# Set output
set terminal gif animate delay 10 loop 0
set output "animation3.gif"


# Configure Axis and title
set size ratio 1
set xrange[0:12.5173]
set yrange[0:12.5173]
set xlabel "X"
set ylabel "Y"
#set title font "Helvetica,14" 

# Setup Color palette
set cbrange [1.5:-0.]
set pal viridis
set tics nomirror

# Loop over frame indices
do for [i=start:end] {
    set title sprintf("Time: %d", i) font "Arial,18"
    plot '../output/velocity.txt' index i u 1:2:3 with image, '../output/velocity.txt' index i u 1:2:($4/(s*sqrt(($4)**2+($5)**2))):($5/(s*sqrt(($4)**2+($5)**2))) every 8:8 with vectors lc -1 filled notitle 
}
