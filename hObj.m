%Defines the handle to allow the optimizations to export additional
%variables, such the as torques, angles and angular velocities of the 7DoF
%upper-limb joints
classdef hObj < handle
   properties
      o=[];
   end
   methods
      function obj=hObj(receivedObject)
         obj.o=receivedObject;
      end
   end
end