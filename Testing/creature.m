classdef creature < handle
    properties
        x
        y
        speed
        fov = 5
        class
        info
    end
    methods
        function moveTo(this, point)
            distance = util.getDistance([this.x this.y], point);
            
            if distance == 0
                return;
            end
            
            directionX = (point(1) - this.x)/distance;
            directionY = (point(2) - this.y)/distance;
            
            if abs(distance)<this.speed
                this.x = point(1);
                this.y = point(2);
                this.speed = 0;
            else
                this.x = this.x+this.speed*directionX;
                this.y = this.y+this.speed*directionY;
            end
        end
        function draw(this, f)
            if this.class == class.predator
                figure(f);
                plot(this.x, this.y, 'r.');
                return;
            end
            if this.class == class.herbivorous
                figure(f);
                plot(this.x, this.y, 'g.');
                return;
            end
        end
        function collectInfo(this)
        end
    end
end