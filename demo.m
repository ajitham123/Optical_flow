% load images:
Im1 = imread('images/shapes_move.png');
Im2 = imread('images/shapes_move_xy135.png');

% determine the individual flow vectors:
graphics = true;
[point, flow_x, flow_y] = get_optical_flow_edges(Im1, Im2, graphics);

% post-process (find FoE, etc.)

