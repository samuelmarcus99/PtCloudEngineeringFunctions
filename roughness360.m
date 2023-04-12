%Function to compute suite of roughness parameters 
function [Rough, maxang, minang, deltaang, directionless, maxrough, minrough, meanprincrough, meanrough, ...
    directionality, netdirectionality, varrough, stdrough, maxdirectionality] = roughness360(data, query, L)
data = data-query; %zero cloud on query point
datat = data'; %transpose raw location data
x = double(data(:,1));
y = double(data(:,2));
z = double(data(:,3));
r = sqrt(sum(data(:,1:2).^2,2)); %find xy distance between each point in the cloud and query point
weight = exp(-r/L);%assign weights to each point
plane = fit([x,y],z,'poly11', 'weight', weight); %fit plane to weighted surface
z = [0, 0, 1];
datat(3,:) = datat(3,:)-plane.p00; %first transformation - remove constant d from f(x,y) = ax+by+d
normal = [plane.p10; plane.p01; -1];
dotted = dot([0; normal(2:3)],z);
denomy = sqrt(0^2+normal(2)^2+normal(3)^2);
firstrot = rotx(acosd(dotted/denomy)); %find angle between plane and y axis
newnormal = firstrot*normal; %rotate plane such that y = 0 - second transformation, to form cf(x,y) = ax
dotted = dot(newnormal,z); 
denomx = sqrt(newnormal(1)^2+newnormal(2)^2+newnormal(3)^2);
beta = acosd(dotted/denomx); %find angle between plane and x axis
secondrot = [cosd(beta) -sind(beta); sind(beta) cosd(beta)]; %reduce dimensionality of rotation - transforms to f(x,y) = 0
newdata = zeros(size(datat));
for num = 1:height(data)
    secondtransform = firstrot*datat(:,num); %apply second transformation to xyz points
    thirdtransform = secondrot*[secondtransform(1); secondtransform(3)]; %apply third transformation to xz points
    newdata(:,num) = [thirdtransform(1); secondtransform(2); thirdtransform(2)]; %reconstruct location matrix - remove the c term from z values
end
sphericaldata=global2localcoord(newdata,"rs"); %convert to spherical coordinates
sphericalangles=sphericaldata(1,:);
sphericalelevation = sphericaldata(3,:).*sind(sphericaldata(2,:)); %reconstruct z values in spherical coordinates
Rough = zeros(1,60); %alter this section based upon the number of bins desired
for i = 1:60
    ang = 3*i - 1.5; %define minimum bin value
    index = ((ang-1.5)<=sphericalangles & sphericalangles<(ang+1.5)|(ang-180-1.5)<=sphericalangles & sphericalangles<(ang-180+1.5));
    heights = sphericalelevation(index);
    Rough(i) = rms(heights); %directional roughness for each point
end
directionless = rms(sphericalelevation); %compute directionless roughness
[maxrough, maxang] = max(Rough); %maximum roughness bin (angle above the horizontal)
[minrough, minang] = min(Rough); %minimum roughness bin (angle above the horizontal)
if abs(maxang-minang) > 30 %sort such that the minimum angle is taken (as opposed to using the clockwise angle)
    deltaang = 3*(60-abs(maxang-minang)); 
else
    deltaang = 3*abs(maxang-minang);
end
meanprincrough = (maxrough+minrough)/2; %mean of the maximum and minimum roughnesses
meanrough = mean(Rough, "omitnan"); %mean of each roughness bin
directionality = Rough./minrough;
maxdirectionality = maxrough/minrough; %directionality of the maximum roughness; from https://www.sciencedirect.com/science/article/pii/S0013795217303551?via%3Dihub#bb0070
netdirectionality = sum(directionality, "omitnan")/60; %sum of the directionality - the greater this value, the greater the anisotropy
varrough = var(Rough, "omitnan"); %variance of the roughness with direction
stdrough = std(Rough, "omitnan"); %standard deviation of the roughness with direction
end