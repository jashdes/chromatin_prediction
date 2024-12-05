function [out] = batchSpotTracking(path, args)

    img = loadImageStack(path);
    out = spotTracking(img, args.template);

end

