clear
close all

% Start
Flt = load('test');
rigidity = 30e9;

% Create stiffness matrices
[stiffk, selfstiff] = stiffnessmatrixk( Flt.cellLength, Flt.cellHeight, ...
    0.5*Flt.cellHeight, rigidity, Flt.nL, Flt.nD );

[stiff, selfstiff1] = stiffnessmatrix( Flt.cellLength, Flt.cellHeight, ...
        0.5*Flt.cellHeight, rigidity, Flt.nL, Flt.nD );

% Calc stress    
stressk = Flt.initStress + slipdeftostressk( Flt.slipDeficit, stiffk );
stress = Flt.initStress + slipdeftostress( Flt.slipDeficit, stiff );
    
figure
imagesc( Flt.slipDeficit ); colorbar;

figure
imagesc( Flt.initStress ); colorbar;

figure
subplot( 3,1,1); imagesc( stress ); colorbar;
subplot( 3,1,2); imagesc( stressk ); colorbar;
subplot( 3,1,3); imagesc( stress - stressk ); colorbar;

[i,j] = find( stressk == min(stressk(:)))
[i,j] = find( stress == min(stress(:)))

figure
imagesc( stress - fliplr(stressk) ); colorbar;
