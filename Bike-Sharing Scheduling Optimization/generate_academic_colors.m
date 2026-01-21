function colors = generate_academic_colors(n)
    if n <= 7
        base_colors = [
            0.000, 0.447, 0.741;
            0.850, 0.325, 0.098;
            0.466, 0.674, 0.188;
            0.929, 0.694, 0.125;
            0.494, 0.184, 0.556;
            0.301, 0.745, 0.933;
            0.635, 0.078, 0.184
        ];
        colors = base_colors(1:n, :);
    else
        colors = zeros(n, 3);
        for i = 1:n
            hue = (i-1) / n;
            colors(i,:) = hsv2rgb([hue, 0.6, 0.8]);
        end
    end
end