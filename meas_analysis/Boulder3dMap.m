function Boulder3dMap()

Glass_wall_BR = [13.199 0 0;13.199 0 7.943;18.368 34.062 7.9430;18.368 34.062 0;13.199 0 0];
Conc_wall_BR = [18.368  34.062 0;18.368  34.062 7.943;0 34.062 7.943;0 34.062 0;18.368  34.062 0];
yw = 33.659;
yt = 34.062;
TX2 = [1.953 1.143 3.364];
TX3 = [10.6440 -1.446 6.2330];


SR1 = [3.046 31.709];
SR2 = [10.035 yw-4.659-1.581];
SR3 = [16.297 yt-12.483];



pos1 = [3.230 2.398 1.917];
pos2 = [3.305 5.668 1.917];
pos3 = [3.641 9.281 1.917];
pos4 = [8.953 yw-1.667-22.229 1.917];
pos5 = [12.204 yw-4.927-18.384 1.917];
pos6 = [12.440 yw-0.772-18.384 1.917];
pos7 = [3.477 14.504 1.917];
pos8 = [3.610 19.037 1.917];
pos9 = [3.349 22.326 1.917];
pos10 = [3.050 25.438 1.917];
pos11 = [2.923 28.601 1.917];

pos12 = [5.076 yw-0.688 1.917];
pos13 = [9.695 yt-0.729 1.917];
pos14 = [14.719 yt-1.646 1.917];
pos15 = [14.955 yt-5.332 1.917];
pos16 = [9.469 yw-4.659-0.918 1.917];

pos17 = [4.572 26.522 1.917];
pos40 = [9.339 26.546 1.917];


pos18 = [9.339 yw-4.659-2.454 1.917];
pos19 = [14.755 yt-8.056 1.917];
pos20 = [14.816 yt-12.208 1.917];
pos21 = [10.161 yw-11.36-1.215 1.917];
pos22 = [4.242 20.992 1.917];
pos23 = [2.276 17.890 1.917];
pos24 = [4.194 14.445 1.917];
pos25 = [9.578 yw-1.141-18.384 1.917];
pos26 = [12.079 yw-1.711-18.384 1.917];
pos27 = [11.331 yw-4.762-18.384 1.917];
pos28 = [8.886 yw-1.008-22.229 1.917];
pos29 = [4.755 10.272 1.917];
pos30 = [4.398 7.635 1.917];
pos31 = [4.815 5.225 1.917];
pos32 = [9.065 yw-2.030-26.510 1.917];
pos33 = [12.181 yw-10.428-18.384 1.917];
pos34 = [11.871 1.678 1.917];
pos35 = [8.784 1.435 1.917];
pos36 = [5.581 1.308 1.917];



antennas = ['SR1';'SR2';'SR3'];
antennax = [SR1(1) SR2(1) SR3(1)];
antennay = [SR1(2) SR2(2) SR3(2)];
labels = ['pos01';'pos02';'pos03';'pos04';'pos05';'pos06';'pos07';'pos08';'pos09';'pos10';'pos11';'pos12';'pos13';'pos14';'pos15';'pos16';'pos17';'pos18';...
            'pos19';'pos20';'pos21';'pos22';'pos23';'pos24';'pos25';'pos26';'pos27';'pos28';'pos29';'pos30';'pos31';'pos32';'pos33';'pos34';'pos35';'pos36'];
% xcoor = [pos1(1) pos2(1) pos3(1) pos4(1) pos5(1) pos6(1) pos7(1) pos8(1) pos9(1) pos10(1) pos11(1) pos12(1) pos13(1) pos14(1) pos15(1) pos16(1)...
%             pos17(1) pos18(1) pos19(1) pos20(1) pos21(1) pos22(1) pos23(1) pos24(1) pos25(1) pos26(1) pos27(1) pos28(1) pos29(1) pos30(1) pos31(1) pos32(1) pos33(1) pos34(1)...
%             pos35(1) pos36(1) pos40(1)];
% ycoor = [pos1(2) pos2(2) pos3(2) pos4(2) pos5(2) pos6(2) pos7(2) pos8(2) pos9(2) pos10(2) pos11(2) pos12(2) pos13(2) pos14(2) pos15(2) pos16(2)...
%             pos17(2) pos18(2) pos19(2) pos20(2) pos21(2) pos22(2) pos23(2) pos24(2) pos25(2) pos26(2) pos27(2) pos28(2) pos29(2) pos30(2) pos31(2) pos32(2) pos33(2) pos34(2)...
%             pos35(2) pos36(2) pos40(2)];
% zcoor = [pos1(3) pos2(3) pos3(3) pos4(3) pos5(3) pos6(3) pos7(3) pos8(3) pos9(3) pos10(3) pos11(3) pos12(3) pos13(3) pos14(3) pos15(3) pos16(3)...
%             pos17(3) pos18(3) pos19(3) pos20(3) pos21(3) pos22(3) pos23(3) pos24(3) pos25(3) pos26(3) pos27(3) pos28(3) pos29(3) pos30(3) pos31(3) pos32(3) pos33(3) pos34(3)...
%             pos35(3) pos36(3) pos40(3)];


% xcoor = [pos1(1) pos2(1) pos3(1) pos4(1) pos5(1) pos6(1) pos7(1) pos8(1) pos9(1) pos10(1) pos11(1) pos12(1) pos13(1) pos14(1) pos15(1) pos16(1)...
%             pos17(1) pos40(1) pos18(1) pos19(1) pos20(1) pos21(1) pos22(1) pos23(1) pos24(1) pos25(1) pos26(1) pos27(1) pos28(1) pos29(1) pos30(1) pos31(1) pos32(1) pos33(1) pos34(1)...
%             pos35(1) pos36(1) ];
% ycoor = [pos1(2) pos2(2) pos3(2) pos4(2) pos5(2) pos6(2) pos7(2) pos8(2) pos9(2) pos10(2) pos11(2) pos12(2) pos13(2) pos14(2) pos15(2) pos16(2)...
%             pos17(2) pos40(2) pos18(2) pos19(2) pos20(2) pos21(2) pos22(2) pos23(2) pos24(2) pos25(2) pos26(2) pos27(2) pos28(2) pos29(2) pos30(2) pos31(2) pos32(2) pos33(2) pos34(2)...
%             pos35(2) pos36(2) ];
% zcoor = [pos1(3) pos2(3) pos3(3) pos4(3) pos5(3) pos6(3) pos7(3) pos8(3) pos9(3) pos10(3) pos11(3) pos12(3) pos13(3) pos14(3) pos15(3) pos16(3)...
%             pos17(3) pos40(3) pos18(3) pos19(3) pos20(3) pos21(3) pos22(3) pos23(3) pos24(3) pos25(3) pos26(3) pos27(3) pos28(3) pos29(3) pos30(3) pos31(3) pos32(3) pos33(3) pos34(3)...
%             pos35(3) pos36(3) ];

xcoor = [pos1(1) pos2(1) pos3(1) pos4(1) pos5(1) pos6(1) pos7(1) pos8(1) pos9(1) pos10(1) pos11(1) pos12(1) pos13(1) pos14(1) pos15(1) pos16(1)...
            pos17(1)  pos18(1) pos19(1) pos20(1) pos21(1) pos22(1) pos23(1) pos24(1) pos25(1) pos26(1) pos27(1) pos28(1) pos29(1) pos30(1) pos31(1) pos32(1) pos33(1) pos34(1)...
            pos35(1) pos36(1) ];
ycoor = [pos1(2) pos2(2) pos3(2) pos4(2) pos5(2) pos6(2) pos7(2) pos8(2) pos9(2) pos10(2) pos11(2) pos12(2) pos13(2) pos14(2) pos15(2) pos16(2)...
            pos17(2)  pos18(2) pos19(2) pos20(2) pos21(2) pos22(2) pos23(2) pos24(2) pos25(2) pos26(2) pos27(2) pos28(2) pos29(2) pos30(2) pos31(2) pos32(2) pos33(2) pos34(2)...
            pos35(2) pos36(2) ];
zcoor = [pos1(3) pos2(3) pos3(3) pos4(3) pos5(3) pos6(3) pos7(3) pos8(3) pos9(3) pos10(3) pos11(3) pos12(3) pos13(3) pos14(3) pos15(3) pos16(3)...
            pos17(3)  pos18(3) pos19(3) pos20(3) pos21(3) pos22(3) pos23(3) pos24(3) pos25(3) pos26(3) pos27(3) pos28(3) pos29(3) pos30(3) pos31(3) pos32(3) pos33(3) pos34(3)...
            pos35(3) pos36(3) ];




        
coord_A1 = [6.971 29.215 3.793];
coord_A2 = [6.971 32.149 3.793];
coord_A3 = [13.665 32.149 3.793];
coord_A4 = [13.665 29.215 3.793];
coord_A5 = [6.971 29.215 0];
coord_A6 = [6.971 32.149 0];
coord_A7 = [13.665 32.149 0];
coord_A8 = [13.665 29.215 0];
obst_A = [coord_A1;coord_A2;coord_A3;coord_A4;coord_A8;coord_A7;coord_A3;coord_A7;coord_A6;coord_A2;coord_A6;coord_A5;coord_A8;coord_A4;coord_A1;coord_A5];
coord_B1 = [6.941 22.530 3.793];
coord_B2 = [6.941 25.481 3.793];
coord_B3 = [13.623 25.481 3.793];
coord_B4 = [13.623 22.530 3.793];
coord_B5 = [6.941 22.530 0];
coord_B6 = [6.941 25.481 0];
coord_B7 = [13.623 25.481 0];
coord_B8 = [13.623 22.530 0];
obst_B = [coord_B1;coord_B2;coord_B3;coord_B4;coord_B8;coord_B7;coord_B3;coord_B7;coord_B6;coord_B2;coord_B6;coord_B5;coord_B8;coord_B4;coord_B1;coord_B5];
coord_C1 = [6.983 15.534 3.807];
coord_C2 = [6.983 18.440 3.807];
coord_C3 = [13.662 18.440 3.807];
coord_C4 = [13.662 15.534 3.807];
coord_C5 = [6.983 15.534 0];
coord_C6 = [6.983 18.440 0];
coord_C7 = [13.662 18.440 0];
coord_C8 = [13.662 15.534 0];
obst_C = [coord_C1;coord_C2;coord_C3;coord_C4;coord_C8;coord_C7;coord_C3;coord_C7;coord_C6;coord_C2;coord_C6;coord_C5;coord_C8;coord_C4;coord_C1;coord_C5];
coord_D1 = [5.937 11.329 4.620];
coord_D2 = [5.937 13.437 4.620];
coord_D3 = [11.403 13.437 4.620];
coord_D4 = [11.403 11.329 4.620];
coord_D5 = [5.937 11.329 0];
coord_D6 = [5.937 13.437 0];
coord_D7 = [11.403 13.437 0];
coord_D8 = [11.403 11.329 0];
obst_D = [coord_D1;coord_D2;coord_D3;coord_D4;coord_D8;coord_D7;coord_D3;coord_D7;coord_D6;coord_D2;coord_D6;coord_D5;coord_D8;coord_D4;coord_D1;coord_D5];
coord_E1 = [7.099 7.300 2.783];
coord_E2 = [7.099 8.927 2.783];
coord_E3 = [10.919 8.927 2.783];
coord_E4 = [10.919 7.300 2.783];
coord_E5 = [7.099 7.300 0];
coord_E6 = [7.099 8.927 0];
coord_E7 = [10.919 8.927 0];
coord_E8 = [10.919 7.300 0];
obst_E = [coord_E1;coord_E2;coord_E3;coord_E4;coord_E8;coord_E7;coord_E3;coord_E7;coord_E6;coord_E2;coord_E6;coord_E5;coord_E8;coord_E4;coord_E1;coord_E5];
coord_F1 = [6.914 2.199 2.789];
coord_F2 = [6.914 4.215 2.789];
coord_F3 = [10.734 4.215 2.789];
coord_F4 = [10.734 2.199 2.789];
coord_F5 = [6.914 2.199 0];
coord_F6 = [6.914 4.215 0];
coord_F7 = [10.734 4.215 0];
coord_F8 = [10.734 2.199 0];
obst_F = [coord_F1;coord_F2;coord_F3;coord_F4;coord_F8;coord_F7;coord_F3;coord_F7;coord_F6;coord_F2;coord_F6;coord_F5;coord_F8;coord_F4;coord_F1;coord_F5];


figure; 
z = ones(size(xcoor)).*1.917;
col = 1:length(xcoor);
plot3(Glass_wall_BR(:,1),Glass_wall_BR(:,2),Glass_wall_BR(:,3),'b')
hold on
plot3(Conc_wall_BR(:,1),Conc_wall_BR(:,2),Conc_wall_BR(:,3),'b')
plot3(obst_A(:,1),obst_A(:,2),obst_A(:,3),'r')
plot3(obst_B(:,1),obst_B(:,2),obst_B(:,3),'r')
plot3(obst_C(:,1),obst_C(:,2),obst_C(:,3),'r')
plot3(obst_D(:,1),obst_D(:,2),obst_D(:,3),'r')
plot3(obst_E(:,1),obst_E(:,2),obst_E(:,3),'r')
plot3(obst_F(:,1),obst_F(:,2),obst_F(:,3),'r')

set(gca,'xaxislocation','top','yaxislocation','right','xdir','reverse')
view(45,75)


% % % surface([xcoor;xcoor],[ycoor;ycoor],[z;z],[col;col],...
% % %         'facecol','no',...
% % %         'edgecol','interp',...
% % %         'linew',2);    
% % % plot3(xcoor,ycoor,zcoor,'k*')

%SO CAN MODIFY THIS TO GET THE CORRECT RESULT...

xcoor2 = xcoor; ycoor2 = ycoor; z2 = z; col2 = col ;

indexvec = col;
%remove these values...
% removervec = [1 2 3];
%removervec = [31 32 33 34 35 36]; %So maybe include 31 as well. ?

%removervec =[19 18  20 21 22 23     removervec ];


removervec = [31 32 33 34 35 36];
removervec = [18 19 20 21 removervec];
removervec = [23 24 25 removervec];




%removervec = [ 21 22 23 24 25 26    removervec]; %I THINK


%removervec = [];


%Now take it over the disjoint set.
indexvec = setdiff(indexvec, removervec);

col2 = col2(indexvec);


%Now adder vec... %we go from 17 - > 22 - > 20 -> 6 + change -> 6
%addvec = [1 20 23 25 50];

%addvec = [];
% % % col3 = col2;
% % % if length(addvec) > 0
% % % for s = 1:length(addvec)
% % %     


% % % %col2 =
% % % col2 = [col2(1:a2-1) addvec(s) col2(a2:end) ]  ;
% % % 
% % % end
% % % end

% 
% if length(removervec) == 0
%Position 6 + change, to avoid an awkward looking line.
pos6change = pos6;
pos6change(1) = pos6change(1) + 2;
xcoor2 = [xcoor2 pos6change(1)];
ycoor2 = [ycoor2 pos6change(2)];
z2 = [z2 pos6change(3)];


addvec = [20];
 [~, a2] = min( abs(  addvec - col2 )  ) ;
 addvecreal = [20 37];
col2 =[col2(1:a2-0) addvecreal col2(a2+1:end) ]  ;



xcoor2 = xcoor2(col2);
ycoor2 = ycoor2(col2);
z2 = z2(col2);



surface([xcoor2;xcoor2],[ycoor2;ycoor2],[z2;z2],[col2;col2],...
        'facecol','no',...
        'edgecol','interp',...
        'linew',2);    
plot3(xcoor2(1:19),ycoor2(1:19),z2(1:19),'k*')
plot3(xcoor2(21:end),ycoor2(21:end),z2(21:end),'k*')


cloudvector = [6  8  13  18  19  29  33];
dx = 0.0; dy = 0.0;
%Run over the cloud locations.
for s = 1:length(cloudvector)
    
    cloudstring = ['Cloud: ', num2str(cloudvector(s))];
    
    text(xcoor(cloudvector(s)) + dx, ycoor(cloudvector(s)) + dy, cloudstring);
    %So we don't double star, but DO star all clouds.
    if ( cloudvector(s) == 29 ||  cloudvector(s) == 33)
    plot3(xcoor(cloudvector(s)) + dx, ycoor(cloudvector(s)) + dy, z(cloudvector(s)), 'k*');
    end
    
end
        
% % % for s = 1:length(xcoor2)
% % %     
% % %     
% % %     
% % %         qtext = num2str(s);
% % %         qtext  = ['Pos ' qtext];
% % %         
% % %         
% % %         
% % %     dx = 0.5 + s/40; dy = 0.1 - s/40;
% % %     text(xcoor2(s)+dx, ycoor2(s)+dy, qtext);
% % %     
% % %     
% % % end




scatter3(TX2(1), TX2(2), TX2(3),'gx')
c = cellstr('TX 2');
dx = 0.5; dy = 0.1;
text(TX2(1)+dx, TX2(2)+dy, TX2(3), c);
scatter3(TX3(1), TX3(2), TX3(3),'gx')
c = cellstr('TX 3');
dx = 0.5; dy = 0.1;
text(TX3(1)+dx, TX3(2)+dy, TX3(3), c);


end

% % % plot(antennax, antennay, 'rx')
% % % for i = 1:length(antennax)
% % %     c = cellstr(antennas);
% % %     dx = 0.5; dy = 0.1;
% % %     text(antennax(i)+dx, antennay(i)+dy, c(i));
% % % end



% % % 
% % % 
% % % for i = 1:length(labels)
% % %     c = cellstr(labels);
% % %     dx = 0.5; dy = 0.1;
% % %     
% % %     
% % %     text(antennax(i)+dx, antennay(i)+dy, c(i));
% % % end
% % % 


