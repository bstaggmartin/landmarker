##Basic Geometric Functions##

#Function to calc distance between two points...
dist.p<-function(p1,p2){
  #Subtract one point from another to get x and y distances, square these elements, sum
  #them, then take square root of the quantity
  return(sqrt(sum((as.numeric(p1)-as.numeric(p2))^2)))
}
#Function to calc circle area with diameter equal to distance between two points...
circ.area<-function(p1,p2){
  #pi*r^2, with r being equal to distance between two points divided by half
  return(pi*((dist.p(p1,p2))/2)^2)
}
#Function to clac area of the triangle formed by three points...
tri.area<-function(p1,p2,p3){
  #Translate vectors going from p1 to p2 and p3 such that they start at 0
  vec1<-p2-p1
  vec2<-p3-p1
  #Use dot product formula to calculate angle between the two vectors
  angle<-acos(sum(vec1*vec2)/(sqrt(sum(vec1^2))*sqrt(sum(vec2^2))))
  #Use base*height/2, with base being equal to distance between p1 and p2, and height being
  #equal to the sin of the angle between the two above vectors multiplied by the distance
  #between p3 and p1 (i.e., the hypotenuse)
  return(dist.p(p1,p2)*dist.p(p1,p3)*sin(angle)/2)
}
#Function used to find angle to transform two points such that they lie in a straight line
#going along the y-axis (with p2 lying in the positive y direction relative to p1)
trans_angle<-function(p1,p2){
  #Find distance between two points, and use that to make a "transformed" point, p.new,
  #that lies the same distance away from p1, but lies away from p1 in the positive y 
  #direction only
  distance<-dist.p(p1,p2)
  p.new<-c(p1[1],p1[2]+distance)
  #Translate vectors going from p1 to p2 and p.new such that they start at 0
  vec1<-p2-p1
  vec2<-p.new-p1
  #Use dot product formula to calculate angle between the two vectors
  #Note: dot product formula calculates magnitude of angular difference, but not direction
  #Use if statement to see if x-coord of p.new is less than that of p2. If so, points must
  #be rotated in clockwise direction, and angular magnitude should remain positive (as is
  #default); otherwise, make angular difference negative to account for transformed
  #position being counter-clockwise relative to original position.
  if(p.new[1]<p2[1]){
    return(acos(sum(vec1*vec2)/(sqrt(sum(vec1^2))*sqrt(sum(vec2^2)))))
  }else{
    return((acos(sum(vec1*vec2)/(sqrt(sum(vec1^2))*sqrt(sum(vec2^2)))))*-1)
  }
}
#Function to find the vector projection of one vector (p1-p0) onto another (p2-p0)
proj<-function(p0,p1,p2){
  p1.new<-p1-p0
  p2.new<-p2-p0
  scalar<-sum(p1.new*p2.new/sqrt(sum(p2.new^2)))
  return(p2.new/sqrt(sum(p2.new^2))*scalar+p0)
}
#Function to rotate a point (p) about another point (p0--the origin by default) by a given angle (ang)
rotate<-function(p,ang,p0=c(0,0)){
  p.new<-p-p0
  return(c(p.new[1]*cos(ang)-p.new[2]*sin(ang),p.new[1]*sin(ang)+p.new[2]*cos(ang))+p0)
}
#Function used to find overlapping area between two circles centered at p1 and p2 and
#radii r1 and r2, respectively
circ.inter<-function(p1,p2,r1,r2,viz=F){
  if(viz==T){
    plot(c(p1[2],p2[2])~c(p1[1],p2[1]),
         xlim=c(min(c(p1[1]-r1,p2[1]-r2)),max(c(p1[1]+r1,p2[1]+r2))),
         ylim=c(min(c(p1[2]-r1,p2[2]-r2)),max(c(p1[2]+r1,p2[2]+r2))),
         xlab="x coord",ylab="y coord",col=c("red","blue"))
    x1<-seq(p1[1]-r1,p1[1]+r1,length.out=100)
    y1<-sqrt(r1^2-(x1-p1[1])^2)
    x2<-seq(p2[1]-r2,p2[1]+r2,length.out=100)
    y2<-sqrt(r2^2-(x2-p2[1])^2)
    lines(p1[2]+y1~x1,col="red");lines(p1[2]-y1~x1,col="red")
    lines(p2[2]+y2~x2,col="blue");lines(p2[2]-y2~x2,col="blue")
  }
  d<-dist.p(p1,p2)
  p1.new<-c(0,0)
  p2.new<-c(d,0)
  c1.l<- -r1
  c1.r<-r1
  c2.l<-d-r2
  c2.r<-d+r2
  if(c1.r==c2.l){
    return("The circles intersect at a single point")
  }else if(c1.r>c2.l & c1.r<c2.r & c2.l>c1.l){
    x.inter<-(d^2+r1^2-r2^2)/(2*d)
    y.inter<-sqrt(r1^2-x.inter^2)
    pos.inter<-c(x.inter,y.inter)
    neg.inter<-c(x.inter,-1*y.inter)
    angle1<-acos(sum(pos.inter*neg.inter)/(sqrt(sum(pos.inter^2))*sqrt(sum(neg.inter^2))))
    angle2<-acos(sum((pos.inter-p2.new)*(neg.inter-p2.new))/
                   (sqrt(sum((pos.inter-p2.new)^2))*sqrt(sum((neg.inter-p2.new)^2))))
    area.contrib1<-r1^2*angle1/2-tri.area(p1.new,pos.inter,neg.inter)
    area.contrib2<-r2^2*angle2/2-tri.area(p2.new,pos.inter,neg.inter)
    return(area.contrib1+area.contrib2)
  }else if(c1.l>=c2.l & c1.r<=c2.r){
    return(pi*r1^2)
  }else if(c1.l<=c2.l & c1.r>=c2.r){
    return(pi*r2^2)
  }else{
    return(0)
  }
}
