# Berechnung der Strommatrix

# Zustände des 36 D Hilbertraums:
Z = cell2mat(struct2cell(load('Hubb_Zust.txt')));
J = zeros(36,36);
for j=1:36; # Spalten (Multiplikation von J mit |j>)
  M = zeros(4,8); # jede erste Komponente negativer Strom, jede zweite Komponente positiver Strom
  for s=1:4; # Alle Elektronen (Einträge in den Zuständen)
    M(:,s*2-1) = Z(:,j);
    M(:,s*2) = Z(:,j);
    # Sprung gegen Stromrichtung:
    if M(s,s*2-1) == 1; # Randbedingung
      M(s,s*2-1) = 4;
    else;
      M(s,s*2-1) -= 1;
    end;
    # Sprung mit Stromrichtung:
    if M(s,s*2) == 4; # Randbedingung
      M(s,s*2) = 1;
    else;
      M(s,s*2) += 1;
    end;
    #if M[s*2,0] != M[s*2,1] or M[s*2+1,2] != M[s*2,3]: # Pauliverbot
    for p = (s*2-1):(s*2); # Summe über alle erzeugten Zustände für das s-te Elektron (jeweils mit und gegen Stromrichtung)
      vzp = (p-s*2+1)*2-1; # -1 oder +1 abhängig von der Stromrichtung

      R = zeros(4,4); # alle 4 Permutationen

      R(:,1) = M(:,p);
      R(1:2,1) = flip(R(1:2,1)); # Spin Up Flip

      R(:,2) = M(:,p);
      R(3:4,2) = flip(R(3:4,2)); # Spin Down Flip

      R(:,3) = M(:,p);
      R(1:2,3) = flip(R(1:2,3)); # Spin Up Flip
      R(3:4,3) = flip(R(3:4,3)); # Spin Down Flip

      R(:,4) = M(:,p); # kein Flip
      
      for l = 1:4; # Vergleich von allen Permutationen des Zustands J|j> mit allen <i|
        vzl = sign(l-2.5); # -1 oder +1 abhängig von der Anzahl an Permutationen
        for i=1:36; # Zeilen (für die Vergleiche)    
          if R(:,l) == Z(:,i); # Vergleich
            J(i,j) += vzp * vzl * 1;
          end;
        end;
      end;
    end;
  end;
end;

# J ist schief-symmetrisch:
for a = 2:36;
  for b = 1:(a-1);
    J(b,a) = -J(a,b);
  end;
end;

J(19:36,19:36)