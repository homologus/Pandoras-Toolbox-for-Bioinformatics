import java.util.*;

import edu.uci.ics.jung.algorithms.shortestpath.DijkstraShortestPath;
import edu.uci.ics.jung.graph.DirectedSparseGraph;


public class My_DFS
{
	private Map<SeqVertex,Integer> _colors;
	public final static int WHITE = 0;
	public final static int BLACK = 1;
	public final static int GRAY = 2;
	public final static int GREEN = 3;
	protected DirectedSparseGraph<SeqVertex,SimpleEdge> _graph;
	protected Set<SeqVertex> _roots;
	
	protected int _time;
	protected Map<SeqVertex,Number> _discovery;
	protected Map<SeqVertex,Number> _finishing;
	
	DijkstraShortestPath<SeqVertex, SimpleEdge> _dp;

	/**
	 * Constructor for the DFS object
	 * Initialize all vertices to WHITE
	 * @param roots 
	 */
	public My_DFS(DirectedSparseGraph<SeqVertex,SimpleEdge> graph) 
	{ 
		_graph = graph;
		initDFS();
	}

	/**
	 * init all variables.
	 */
	public void initDFS()
	{
		_colors = new HashMap<SeqVertex,Integer>();
		for (SeqVertex v : _graph.getVertices())
			_colors.put(v,WHITE);
		
		_time=0;
		_discovery = new HashMap<SeqVertex,Number>();
		_finishing = new HashMap<SeqVertex,Number>();
		_dp = new DijkstraShortestPath<SeqVertex, SimpleEdge>(_graph);
	}
	
	public void runDFS()
	{
		initDFS();
		for (SeqVertex v : _graph.getVertices())
			if (_graph.inDegree(v)==0 && getColor(v)==WHITE) {
				visitVertex(v);
			}
		
		//update_start_node_depths_bottom_up();
		
//		updateCircularity();
	}
	
	
	public void update_start_node_depths_bottom_up() {
		
		for (SeqVertex v : _graph.getVertices()) {
			if (_graph.inDegree(v) == 0 && _graph.outDegree(v) == 1) {
				// start node that is itself unbranched.
				// walk it until we find a branch point.
				
				List<SeqVertex> node_list = new ArrayList<SeqVertex>();
				node_list.add(v);
				
				SeqVertex next_v = v;
				
				while (_graph.outDegree(next_v) == 1 && _graph.inDegree(next_v) == 1) {
					node_list.add(next_v);
					next_v = _graph.getSuccessors(next_v).iterator().next();
				}
				
				if (_graph.getSuccessors(next_v).iterator().hasNext()) { // next is a branched (before or after) node.
					int branched_node_depth = _graph.getSuccessors(next_v).iterator().next()._node_depth;
					Collections.reverse(node_list);
					for (SeqVertex u : node_list) {
						branched_node_depth--;
						u._node_depth = branched_node_depth;
					}
				}
				
			}
		}
		
	}
	
	
	
	
	
//	private void updateCircularity() {
//		DijkstraShortestPath<SeqVertex, SimpleEdge> dp = new DijkstraShortestPath<SeqVertex, SimpleEdge>(_graph);
//
//		for (SeqVertex v : _graph.getVertices())
//			if (v.isInCircle())
//			{
//				assert(dp.getDistance(v, v)!=null); //v is reachable from itself
//				for (SimpleEdge e : dp.getPath(v, v))
//				{
//					e.setInCircle(true);
//					_graph.getDest(e).setInCircle(true);
//				}
//			}
//	}

	/**
	 * Visit a vertex:
	 * color it gray, and explore all its descendants
	 */
	private void visitVertex(SeqVertex v)
	{
		setColor(v, GRAY);
		_time++;
		// mark the discovery time of this vertex
		_discovery.put(v, _time);
		v.setDFS_DiscoveryTime(_time);
		
		// determine node depth based on min (parent node depth) + 1
		Integer node_depth = null;
		for (SeqVertex p : _graph.getPredecessors(v)) {
			if (node_depth == null)
				node_depth = p._node_depth;
			else 
				node_depth = Math.max(node_depth, p._node_depth);
				
		}
		if (node_depth == null)
			node_depth = 0;
		else
			node_depth += 1;
		
		v._node_depth = node_depth;
		
		for (SeqVertex u : _graph.getSuccessors(v))
		{
			if (getColor(u)==WHITE)
				visitVertex(u);
			else if (getColor(u)==GRAY){ // we have reached a circle.
				u.setInCircle(true);
				assert(_dp.getDistance(u, v)!=null); //v is reachable from itself
				if (v.equals(u)) //self loop
					_graph.findEdge(v, v).setInCircle(true);
				else //other loops
					for (SimpleEdge e : _dp.getPath(u, v))
					{
						e.setInCircle(true);
						_graph.getDest(e).setInCircle(true);
					}

			}
		}
		
		setColor(v, BLACK);
		_time++;
		// mark the finishing time of this vertex
		_finishing.put(v, _time);
		v.setDFS_FinishingTime(_time);
	}
	
	/**
	 * Returns the color of the given vertex	
	 */
	public int getColor(SeqVertex v)
	{
		return  _colors.get(v);
	}

	/**
	 * Sets the color of the given vertex	
	 */
	public void setColor(SeqVertex v, Integer color)
	{
		_colors.put(v,color);
	}

	/**
	 * return the finishing time for each vertex
	 * @return
	 */
	public Map<SeqVertex, Number> getFinishing() {
		return _finishing;
	}



}



