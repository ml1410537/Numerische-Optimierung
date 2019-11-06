param pi := 3.14159265358979;
# Schiffsparameter
param c1 := -0.26;
param c2 := 0.2; 
param c3 := -1.87;
param c4 := 0.6;
# Geschwindigkeit
param v := 3.5;
# Stellgrößenbeschränkungen
param umin := -15 * pi/180;
param umax := 15 * pi/180;
# Endzeit
param tf := 20;

set Nx := {1..4};
# Anfangsbedingungen    
param x0 {i in Nx};
let x0[1] := 0;
let x0[2] := 0;
let x0[3] := 0;
let x0[4] := 0;
# (Gewünschte) Endbedingungen
param xf {i in Nx};
let xf[1] := 0;
let xf[2] := 0;
let xf[3] := 0;
let xf[4] := 10;

# Gewichtungen des Endzustands      
param S {i in Nx};
let  S[1] := 1;
let  S[2] := 1;
let  S[3] := 1;
let  S[4] := 1;
# Gewichtungen der Zustände im Integralanteil
param Q {i in Nx};
let Q[1] := 1;
let Q[2] := 1;
let Q[3] := 1;
let Q[4] := 1;
# Gewichtung der Stellgröße im Integralanteil
param R := 0.1;

# Anzahl Diskretisierungspunkte
param N   := 59;
# Konstante Schrittweite h = t[k+1] - t[k]
param h   := tf/N;
set Nt  := {0..N};
set Nt1 := {0..N-1};

# Optimierungsvariablen
# ---------------------
# Stellgrößen mit Beschränkungen
var u{Nt} <= umax, >= umin;
# Zustände
var x{Nx, Nt};

# Diskretisiertes Kostenfunktional (Trapezregel)
# ----------------------------------------------
minimize cost:
    # Gewichtung des Endzustands
    sum{i in Nx}( ??? ) +
    # Gewichtung des Integralanteils
    h/2 * sum{k in Nt1}(
        sum{i in Nx}( ??? )
    );

# Diskretisierte Systemgleichungen (Trapezregel)
# ----------------------------------------------
subject to x1_diskr {k in Nt1} : 
    x[1,k+1] - x[1,k] = h/2 * ( ??? );

subject to x2_diskr {k in Nt1} :  
    x[2,k+1] - x[2,k] = h/2 * ( ??? );

subject to x3_diskr {k in Nt1} : 
    x[3,k+1] - x[3,k] = h/2 * ( ??? );

subject to x4_diskr {k in Nt1} : 
    x[4,k+1] - x[4,k] = h/2 * ( ??? );

# Anfangsbedingungen
# ------------------
subject to x_init {i in Nx} : x[i,0] = x0[i] ;

# Loesung und Ausgabe
# -------------------
option solver '/path/to/solver';
solve;
display cost ;

# Ausgabe in Datei
printf {k in Nt}: "%10.5f %10.5f %10.5f %10.5f %10.5f %10.5f\n", 
                  k*h, x[1,k], x[2,k], x[3,k], x[4,k], u[k] > schiff_AMPL_sol.txt;
printf "Loesung wurde in 'schiff_AMPL_sol.txt' gespeichert.\n\n";
