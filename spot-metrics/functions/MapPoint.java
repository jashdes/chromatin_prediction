public class MapPoint implements Comparable<MapPoint> {

	// Wrapper class to contain point indeces and their distances
	// for use in a PriorityQueue

	public int index;
	public double value;

	public MapPoint(int i, double v){
		index = i;
		value = v;
	}

	@Override
	public int compareTo(MapPoint o){
		if(this.value == o.value){
			return 0;
		} else {
			return this.value > o.value ? 1 : -1;
		}
	}

}
