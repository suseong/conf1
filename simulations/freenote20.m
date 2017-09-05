for kkkk=1:length(segInfo)
    for kkk=1:length(segInfo{kkkk})
        funnelIdx = segInfo{kkkk}{kkk}{2};
        pos{kkkk}{kkk} = segInfo{kkkk}{kkk}{1};
        if funnelIdx(1) == 1
            chad = toEllipsoid(funnel_s{funnelIdx(2)}(:,funnelIdx(3)));
        else
            chad = toEllipsoid(funnel_b{funnelIdx(2)}(:,funnelIdx(3)));
            [fx,fy,fz] = sphere(10);
            fx = fx/sqrt(Px(1,1));
            fy = fy/sqrt(Px(2,2));
            fz = fz/sqrt(Px(3,3));
            Px
        end
        Pp = chad(1:3,1:3); Ppv = chad(1:3,4:6); Pv = chad(4:6,4:6);
        Px = Pp-Ppv^2*pinv(Pv);
        funnelSize{kkkk}{kkk} = sqrt(diag(Px));
    end
end