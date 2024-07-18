C = [ ...
    204, 204, 204; % Gray
    68 , 119, 170; % Blue
    238, 102, 119; % Red
    204, 187, 68 ; % Yellow
    102, 204, 238; % Cyan
    170, 51 , 119; % Purple
    34 , 136, 51 ; % Green
    ];

C  = C/255;

tolgray = C(1, :); 
tolblue = C(2, :);
tolred = C(3, :);
tolyellow = C(4, :);
tolcyan = C(5, :);
tolpurple = C(6, :);
tolgreen = C(7, :);

save("colors.mat", "tolgray", "tolblue", "tolred", "tolyellow", "tolcyan", "tolpurple","tolgreen"); 