import qupath.lib.objects.PathObjects
import qupath.lib.roi.ROIs
/**
Use this script to create a centroid reference point in the annotations to be used in tissue reconstruction
Use the centroid in combination with a rectangle annotation to aid in placing cardinal reference points that are perpendicular
Parameters to set: 
1. classes used for tissue annotations, change as needed
2. Size of the ellipse to create for tissue centroid in pixels
*/
def classes = [getPathClass("Class1"),getPathClass("Class2"),getPathClass("Class3"),getPathClass("Class4")]
double height = 500
double width = 500
def centroidClass = 'className' //class to use for the centroid object
def name = 'objectName' //name to give the centroid object

//find all annotations matching the classes above
def imageData = getCurrentImageData()
def hierarchy = imageData.getHierarchy()
def tissues = []

classes.each { class ->
   def annos = getAnnotationObjects().findAll {it.getPathClass() == class}
   tissues << annos;
}

def newAnnotations = []
println "Found annotations:"
println tissues
tissues.each { tissue ->
   def roi = tissue.getROI()
   double x = roi.getCentroidX()-250
   double y = roi.getCentroidY()-250
   def object = ROIs.createEllipseROI(x,y,height,width)
   def newAnno = PathObjects.createAnnotationObject(object,getPathClass(centroidClass))
   newAnno.setName(name)
   newAnnotations << newAnno;
}

println 'adding new objects'
addObjects(newAnnotations)