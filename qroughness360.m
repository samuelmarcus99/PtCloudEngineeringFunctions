%Function to determine roughness projected onto a plane, fitted to a local quadric
%Identical to roughness360 to line 33
function [qRough, qmaxang, qminang, qdeltaang, qdirectionless, qmaxrough, qminrough, qmeanprincrough, qmeanrough, ...
    qdirectionality, qnetdirectionality, qvarrough, qstdrough, qmaxdirectionality] = qroughness360(data, query, L)
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
z = sphericaldata(3,:).*sind(sphericaldata(2,:)); %reconstruct z values in spherical coordinates
x = sphericaldata(3,:).*cosd(sphericaldata(1,:)); %reconstruct x and y values in spherical coordinates (to allow for quadric fitting)
y = sphericaldata(3,:).*sind(sphericaldata(1,:));
quad = fit([x',y'], z', 'poly22', 'weight', weight); %fit quadric using reconstructed xy values
qRough = zeros(1,60);
for i = 1:60
    ang = 3*i - 1.5;
    ixx = ((ang-1.5)<=sphericalangles & sphericalangles<(ang+1.5)|(ang-180-1.5)<=sphericalangles & sphericalangles<(ang-180+1.5));
    heights = z(ixx) - feval(quad, x(ixx), y(ixx)); %take heights by evaluating quadric at each point and subtracting z values from quadric
    qRough(i) = rms(heights); 
end
fun = z - feval(quad,x,y); %take heights by evaluating quadric at each point and subtracting z values from quadric
qdirectionless = rms(fun); %from here, once again identical to roughness360
[qmaxrough, qmaxang] = max(qRough); %maximum roughness bin (angle above the horizontal)
[qminrough, qminang] = min(qRough); %minimum roughness bin (angle above the horizontal)
if abs(qmaxang-qminang) > 30
    qdeltaang = 3*abs(60-(qmaxang-qminang));
else
    qdeltaang = 3*abs(qmaxang-qminang);
end
qmeanprincrough = (qmaxrough+qminrough)/2; %mean of the maximum and minimum roughnesses
qmeanrough = mean(qRough, "omitnan"); %mean of each roughness bin
qdirectionality = qRough./qminrough;
qmaxdirectionality = qmaxrough/qminrough; %directionality of the maximum roughness; from https://www.sciencedirect.com/science/article/pii/S0013795217303551?via%3Dihub#bb0070
qnetdirectionality = sum(qdirectionality, "omitnan")/60; %net sum of the directionality - the greater this value, the greater the anisotropy
qvarrough = var(qRough, "omitnan"); %variance of the roughness with direction
qstdrough = std(qRough, "omitnan"); %standard deviation of the roughness with direction
end