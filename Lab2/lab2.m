clear; clc;


%% Mandelbrot Set
x = linspace(-2, 2, 5000);
y = linspace(-2, 2, 5000);

insideCount = 0;
flag = 0;

for i = 1:5000
    for j = 1:5000
        if (x(i)^2 + y(j)^2 <= 4)
            insideCount = insideCount+1;
            Inside_Point(insideCount) = complex(x(i), y(j));
        end
    end
end

plot(Inside_Point, '.', 'MarkerEdgeColor', 'k');

hold on;
Original_Inside_Point = Inside_Point;

for iteration = 1:600
 
    len = insideCount;
    outsideCount = 0;
    outside_Point = [];
    insideCount = 0;
    temp_Inside_Point = [];
    temp_Original_Inside_Point = [];

    for i = 1:len
        
        % Mandelbrot set
        result = Inside_Point(i)^2 + Original_Inside_Point(i);

        % outside
        if(abs(result) > 2)
            outsideCount = outsideCount + 1;
            outside_Point(outsideCount) = Original_Inside_Point(i);
            if(iteration == 600 && flag == 0)
                a = Original_Inside_Point(i);
                flag = 1;
            end
        % inside
        else
            insideCount = insideCount + 1;
            temp_Inside_Point(insideCount) = result;
            temp_Original_Inside_Point(insideCount) = Original_Inside_Point(i);
        end

    end

    if mod(iteration, 3) == 1
        plot(outside_Point, '.', 'MarkerEdgeColor', 'r');
    elseif mod(iteration, 3) == 2
        plot(outside_Point, '.', 'MarkerEdgeColor', 'g');
    else
        plot(outside_Point, '.', 'MarkerEdgeColor', 'b');
    end

    Inside_Point = temp_Inside_Point;
    Original_Inside_Point = temp_Original_Inside_Point;
end

%% print
frame_count = 1;
xx1 = -2; xx2 = 2; yy1 = -2; yy2 = 2; r = 1;
midpoint_x = -1.0862;
midpoint_y = -0.2365;
for iteration = 1:620
    axis([xx1,xx2,yy1,yy2]);
    xx1 = midpoint_x-(midpoint_x-(-2))/r;
    xx2 = midpoint_x-(midpoint_x-2)/r;
    yy1 = midpoint_y-(midpoint_y-(-2))/r;
    yy2 = midpoint_y-(midpoint_y-2)/r;
    r = r+0.04;
    set(gcf,'position',[0.1,0.1, 800, 800]);
    frame_array(frame_count) = getframe;
    frame_count = frame_count + 1;
end

v = VideoWriter('lab2');
v.FrameRate = 10;
open(v);
writeVideo(v, frame_array);
close(v);