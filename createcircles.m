clc;
clear;

rect_width = 20;     
rect_height = 40;     
area_ratio = 0.5;     
rmin=0.08;
rmax=0.12;
part_O=2;
P1_1=part_O-part_O*0.1;
P1_2=part_O+(rect_height-3*part_O)/2+part_O*0.1;
P2_1=((rect_height-3*part_O)/2+2*part_O)-part_O*0.1;
P2_2=(2*(rect_height-3*part_O)/2+2*part_O)+part_O*0.1;
freq_data = [2, 32, 37, 21, 7, 1];
x_midpoints = ([55, 70, 100, 130, 160, 190]) * 0.001; 
mu = mean(log(x_midpoints)); 
sigma = std(log(x_midpoints)); 

circles_part_O = generate_circles_in_area(0, part_O, rect_width, area_ratio, mu, sigma, freq_data, [],rmin,rmax);
circles_part_O1 = circles_part_O;

circles_part_O(:,2) = circles_part_O(:,2) + (rect_height-3*part_O)/2+part_O; 
circles_Part1_O = circles_part_O;

circles_part_O(:,2) = circles_part_O(:,2) + (rect_height-3*part_O)/2+part_O; 
circles_part2_O = circles_part_O;

all_circles = [circles_part_O1; circles_Part1_O; circles_part2_O];

circles_part1 = generate_circles_in_area(P1_1, P1_2, rect_width, area_ratio, mu, sigma, freq_data, all_circles,rmin,rmax);
circles_part2 = generate_circles_in_area(P2_1, P2_2, rect_width, area_ratio, mu, sigma, freq_data, all_circles,rmin,rmax);

all_circles = [all_circles; circles_part1; circles_part2];

% figure(1);
% hold on;
% plot_circles(circles_part_O1, 0, 'r');     
% plot_circles(circles_Part1_O, 0, 'g');        
% plot_circles(circles_part2_O, 0, 'b');        
% plot_circles(circles_part1, 0, 'm');    
% plot_circles(circles_part2, 0, 'c');     
% axis equal;
% xlim([0 rect_width]);
% ylim([0 rect_height]);
% xlabel('Width (µm)');
% ylabel('Height (µm)');
% title('Generated Circles with Copy');
% hold off;

mphopen('G:\learn\WPSstorefile\finalpaper\paper3\crateparticles\DEMO.mph')
circles=all_circles;
W=[num2str(rect_width) '[um]'];
H=[num2str(rect_height) '[um]'];
model.param.set('H', H);
model.param.set('W', W);
model.param.set('W', W);
model.param.set('area_ratio', num2str(area_ratio));
[m n]=size(circles);
H_sep=20;
H_cur=10;
cc=[];
ccc=[];
for ii=1:m
c=['C',num2str(ii)];
model.geom('geom1').create(c, 'Circle');
model.geom('geom1').feature(c).set('r', circles(ii,3));
model.geom('geom1').feature(c).set('pos', [circles(ii,1) circles(ii,2)]);
model.geom('geom1').feature(c).set('selresult', true);
ccc=[ccc ' ' c];
cc=[cc ' ' '''' c ''''];
end
trimmedString = strtrim(cc);
trimmedString1 = strtrim(ccc);
model.geom('geom1').create('sel1', 'ExplicitSelection');
for jj=1:m
c=['C',num2str(jj)];
model.geom('geom1').feature('sel1').selection('selection').set(c, 1);
end 
model.geom('geom1').create('r_1', 'Rectangle');
model.geom('geom1').feature('r_1').set('size', [rect_width rect_height]);
model.geom('geom1').feature('r_1').set('pos', [0 0]);
model.geom('geom1').create('r_2', 'Rectangle');
model.geom('geom1').feature('r_2').set('size', [rect_width H_sep]);
model.geom('geom1').feature('r_2').set('pos', [0 -H_sep]);
model.geom('geom1').create('r_3', 'Rectangle');
model.geom('geom1').feature('r_3').set('size', [rect_width H_cur]);
model.geom('geom1').feature('r_3').set('pos', [0 rect_height]);
model.geom('geom1').create('uni1', 'Union');
cellArray = regexp(trimmedString, 'C\d+', 'match');
cellArray = cellfun(@(x) {x}, cellArray);
model.geom('geom1').feature('uni1').selection('input').set(cellArray);
model.geom('geom1').feature('uni1').set('selresult', true);
model.geom('geom1').create('par1', 'Partition');
model.geom('geom1').feature('par1').selection('input').set({'r_1'});
model.geom('geom1').feature('par1').selection('tool').named('uni1');
model.geom('geom1').feature('par1').set('selresult', true);
model.geom('geom1').create('intsel1', 'IntersectionSelection');
model.geom('geom1').feature('intsel1').set('input', {'uni1'});
model.geom('geom1').feature('intsel1').set('entitydim', 1);
model.geom('geom1').run('fin');
figure(2)
%mphgeom(model)

function circles = generate_circles_in_area(z1, z2, rect_width, area_ratio, mu, sigma, freq_data, existing_circles, rmin, rmax)
       area_rect = rect_width * (z2 - z1);
    area_circles = area_rect * area_ratio;
    max_circles = floor(area_circles / (pi * (0.15)^2));
    
    circles = [];

    while size(circles, 1) < max_circles
               radius = lognrnd(mu, sigma); 
        if radius < rmin || radius > rmax  
            continue;
        end
        x = rand() * (rect_width - 2 * radius) + radius; 
        y = rand() * (z2 - z1 - 2 * radius) + (z1 + radius); 

        if is_valid_circle(circles, [x, y], radius) && is_valid_circle(existing_circles, [x, y], radius)
            circles = [circles; x, y, radius]; 
        end
    end
end



function is_valid = is_valid_circle(existing_circles, new_circle, radius)
     is_valid = true;
    min_distance = 0.1 * radius; 

    for i = 1:size(existing_circles, 1)
        if norm(existing_circles(i, 1:2) - new_circle) < (radius + existing_circles(i, 3)) + min_distance
            is_valid = false;
            return;
        end
    end
end

function plot_circles(circles, z_height, color)
     theta = linspace(0, 2*pi, 100);
    for i = 1:size(circles, 1)
        x_circle = circles(i, 1) + circles(i, 3) * cos(theta);
        y_circle = circles(i, 2) + circles(i, 3) * sin(theta) + z_height;
        fill(x_circle, y_circle, color, 'FaceAlpha', 0.5);
    end
end