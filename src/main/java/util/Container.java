package util;


import java.util.*;


/**
 * Created by 左岗岗 on 2017/1/8.
 */
public class Container<T> {

    public ArrayList<T> creatConstantList(T value, int size){
        ArrayList<T> constantList = new ArrayList<T>();
        for (int i = 0; i < size; i++) {
            constantList.add(value);
        }
        return constantList;
    }

    public ArrayList<T> toArrayList(Set<T> target){
        ArrayList<T> result = new ArrayList<T>();
        for (T s:target){
            result.add(s);
        }
        return result;
    }

    public ArrayList<T> withoutRepeat(ArrayList<T> arrayList){
        Set<T> notRepeat = new TreeSet<T>();
        for (int i = 0; i < arrayList.size(); i++) {
            notRepeat.add(arrayList.get(i));
        }

        ArrayList<T> result = new ArrayList<T>();
        for (T s : notRepeat) {
            result.add(s);
        }
        return result;
    }

    public static double[][] as2DDoubleArray(String[][] arr) {
        double[][] doubArr = new double[arr.length][arr[0].length];
        for (int i = 0; i < arr.length; i++) {
            for (int j = 0; j < arr[i].length; j++) {
                doubArr[i][j] = Double.parseDouble(arr[i][j].trim());
            }
        }
        return doubArr;
    }

    public static double[] asDoubleArray(String[] arr){
        double[] result = new double[arr.length];
        for (int i = 0; i < result.length; i++) {
            result[i]=Double.valueOf(arr[i].trim());
        }
        return result;
    }

    public static double[] toDoubleArray(List<Double> list){
        double[] res = new double[list.size()];
        for (int i = 0; i < res.length; i++) {
            res[i] = list.get(i);
        }
        return res;
    }



    public static ArrayList<String> toStringArrayList(ArrayList<Double> arrayList) {
        ArrayList<String> result = new ArrayList<String>();
        for (int i = 0; i < arrayList.size(); i++) {
            result.add(String.valueOf(arrayList.get(i)).trim());
        }
        return result;
    }

    public static ArrayList<String> toStringArrayList(double[] arr) {
        ArrayList<String> result = new ArrayList<String>();
        for (int i = 0; i < arr.length; i++) {
            result.add(String.valueOf(arr[i]));
        }
        return result;
    }

    public ArrayList<T> toArrayList(T[] array){
        ArrayList<T> arrayList = new ArrayList<T>();
        for (int i = 0; i < array.length; i++) {
            arrayList.add(array[i]);
        }
        return arrayList;
    }
    public static ArrayList<Double> toArrayList(double[] array){
        ArrayList<Double> arrayList = new ArrayList<Double>();
        for (int i = 0; i < array.length; i++) {
            arrayList.add(array[i]);
        }
        return arrayList;
    }
    public static HashMap<Object,Object> upcastToObject(HashMap<Object,Double> hashMap){
        HashMap<Object,Object> newHashMap = new HashMap<Object, Object>();
        for (int i = 0; i < hashMap.keySet().size(); i++) {
            newHashMap.put(hashMap.keySet().toArray()[i],hashMap.get(hashMap.keySet().toArray()[i]));
        }
        return newHashMap;
    }

    public static double[][] indexArray(double[] arr,int startIndex){
        double[][] result=new double[arr.length][2];
        for (int i = 0; i < arr.length; i++) {
            result[i][0]=i+startIndex;//index
            result[i][1]=arr[i];//value
        }
        return result;
    }

    public static double[][] indexAndListCommon(double common,int num,int startIndex){
        double[][] result=new double[num][2];
        for (int i = 0; i < num; i++) {
            result[i][0]=i+startIndex;
            result[i][1]=common;
        }
        return result;
    }

    public static double[][] indexAndListOppositeCommon(double common,int num,int startIndex){
        double[][] result=new double[num][2];
        for (int i = 0; i < num; i++) {
            result[i][0]=i+startIndex;
            result[i][1]=-common;
        }
        return result;
    }








}
