/*
Attempts to create a minimum spanning tree, given a set of points and and
a set of edges. The initial pool of edges is shuffled, and then the function
iterates through the pool of edges, adding edges until either all points
are connected or the function reaches the end of the list of edges.
*/

vector<Edge> random_mst(set<Point> points, vector<Edge> edges){

    unordered_map<Point, int> ufmap;
    vector<Edge> kept;
    vector<int> canonical(points.size(), -1);
    int index = 0, num_roots = points.size(), root1, root2;

    for(Point p : points){
        ufmap.insert(make_pair(p, index++));
    }

    shuffle_edges(edges);

    for(Edge next : edges){
        root1 = find(canonical, ufmap[next.p1]);
        root2 = find(canonical, ufmap[next.p2]);
        if(root1 != root2){
            makeunion(canonical, root1, root2);
            kept.push_back(next);
            if(--num_roots == 1) break;
        }
    }

    return kept;
}
