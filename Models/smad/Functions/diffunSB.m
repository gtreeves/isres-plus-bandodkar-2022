function dydt = diffunSB(t,y,k)


R = y(1);
TGFb = y(2);
Ract = y(3);
Rinact = y(4);
SB  = y(5);
S2c = y(6);
Gc = y(7);
pS2c = y(8);
pGc = y(9);
S4c = y(10);
S24c = y(11);
G4c = y(12);
S22c = y(13);
G2c = y(14);
GGc = y(15);
S2n = y(16);
Gn = y(17);
pS2n = y(18);
pGn = y(19);
S4n = y(20);
S24n = y(21);
G4n = y(22);
S22n = y(23);
G2n = y(24);
GGn = y(25);

PPase = 1;

dR = -k(1)*R*TGFb;
dTGFb = -k(1)*R*TGFb;
dRact = k(1)*R*TGFb-k(3)*Ract*SB+k(10)*Rinact;
dRinact = k(3)*Ract*SB-k(10)*Rinact;
dSB = -k(3)*Ract*SB+k(10) *Rinact;
dS2c = k(8)*S2n-k(7)*S2c-k(2)*S2c*Ract;
dGc = k(8)*Gn-k(7)*Gc-k(2)*Gc*Ract;
dpS2c = k(8)*pS2n-k(7)*pS2c+k(2)*S2c*Ract-k(4)*pS2c*(S4c+2*pS2c+pGc)+k(9)*(S24c+2*S22c+G2c);
dpGc = k(8)*pGn - k(7)*pGc+ k(2)*Gc*Ract-k(4)*pGc*(S4c+pS2c+2*pGc)+k(9)*(G4c+G2c+2*GGc); 
dS4c = k(7)*S4n-k(7)*S4c-k(4)*S4c*(pS2c+pGc)+k(9)*(S24c+G4c);
dS24c = k(4)*S4c*pS2c - k(9)*S24c - k(7)*k(6)*S24c;
dG4c = k(4)*pGc*S4c-k(9)*G4c-k(7)*k(6)*G4c;
dS22c = k(4)*pS2c*pS2c-k(9)*S22c-k(7)*k(6)*S22c;
dG2c =k(4)*pGc*pS2c-k(9)*G2c-k(7)*k(6)*G2c;
dGGc =k(4)*pGc*pGc-k(9)*GGc-k(7)*k(6)*GGc;
dS2n = k(7)*S2c - k(8)*S2n + k(5)*pS2n*PPase;
dGn =k(7)*Gc-k(8)*Gn+k(5)*pGn*PPase;
dpS2n = k(7)*pS2c - k(8)*pS2n - k(5)*pS2n*PPase - k(4)*pS2n*(S4n+2*pS2n+pGn)+k(9)*(S24n+2*S22n+G2n);
dpGn = k(7)*pGc - k(8)*pGn- k(5)*pGn*PPase-k(4)*pGn*(S4n+pS2n+2*pGn)+k(9)*(G4n+G2n+2*GGn); 
dS4n = k(7)*S4c- k(7)*S4n- k(4)*S4n*(pS2n+pGn)+ k(9)*(S24n+G4n);
dS24n = k(4)*S4n*pS2n - k(9)*S24n + k(7)*k(6)*S24c;
dG4n = k(4)*pGn*S4n-k(9)*G4n+k(7)*k(6)*G4c;
dS22n = k(4)*pS2n*pS2n-k(9)*S22n+k(7)*k(6)*S22c;
dG2n = k(4)*pGn*pS2n-k(9)*G2n+k(7)*k(6)*G2c;
dGGn = k(4)*pGn*pGn-k(9)*GGn+k(7)*k(6)*GGc;

dydt = [dR;dTGFb;dRact;dRinact;dSB;dS2c;dGc;dpS2c;dpGc;dS4c;dS24c;dG4c;dS22c;dG2c;dGGc;dS2n;dGn;dpS2n;dpGn;dS4n;dS24n;dG4n;dS22n;dG2n;dGGn];

end