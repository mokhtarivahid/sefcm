
data1 = input( 'Enter the name of cluster file: ', 's' );
data2 = input( 'Enter the name of centeroids file: ', 's' );

fcmdata = load( data1 );
center = load( data2 );

color = {'m' 'c' 'r' 'g' 'b' 'k' 'y' 'w'};

hold on
label = 0;
for i = 1 : max(fcmdata(:,3))+1
    index = find( fcmdata(:,3) == label );
    line( fcmdata(index, 1), fcmdata(index, 2), 'linestyle', 'none', 'marker', 'o', 'color', color{i} );
    plot( center(i,1), center(i,2), 'ko','markersize',12,'LineWidth', 2 )
    label = label + 1;
end
