function u = generate_u(u0, interval, vals)
	len = (size(vals, 2)+1)*interval;
    u = u0*ones(1, len);
    for i = 1:size(vals, 2)
		u(i*interval:end) = vals(i);
    end
end

