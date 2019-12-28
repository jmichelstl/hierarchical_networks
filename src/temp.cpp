bool intersection(Edge e1, Edge e2){
    double slope1, slope2, x;
    Point p1 = e1.p1, p2 = e1.p2, p3 = e2.p1, p4 = e2.p2;
    bool inxrange1, inxrange2;

    slope1 = p1.x == p2.x ? BIG_SLOPE : (p2.y - p1.y)/(p2.x - p1.x);
    slope2 = p3.x == p4.x ? BIG_SLOPE : (p4.y - p3.y)/(p4.x - p3.x);

    x = (p1.y - p3.y + slope2*p3.x - slope1*p1.x)/(slope2 - slope1);

    inxrange1 = x >= min(p1.x,p2.x) && x <= max(p1.x, p2.x);
    inxrange2 = x >= min(p3.x,p4.x) && x <= max(p3.x, p4.x);

    return inxrange1 && inxrange2;
}
