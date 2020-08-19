use super::compare_segments::compare_segments;
use super::compute_fields::compute_fields;
use super::helper::{BoundingBox, Float};
use super::possible_intersection::possible_intersection;
use super::sweep_event::SweepEvent;
use super::Operation;
use array_stump::ArrayStump;
use std::collections::BinaryHeap;
use std::rc::Rc;

#[cfg(feature = "debug-booleanop")]
use super::sweep_event::JsonDebug;

pub fn subdivide<F>(
    event_queue: &mut BinaryHeap<Rc<SweepEvent<F>>>,
    sbbox: &BoundingBox<F>,
    cbbox: &BoundingBox<F>,
    operation: Operation,
) -> Vec<Rc<SweepEvent<F>>>
where
    F: Float,
{
    let mut sweep_line = ArrayStump::<Rc<SweepEvent<F>>, _>::new(compare_segments);
    let mut sorted_events: Vec<Rc<SweepEvent<F>>> = Vec::new();
    let rightbound = sbbox.max.x.min(cbbox.max.x);

    while let Some(event) = event_queue.pop() {
        #[cfg(feature = "debug-booleanop")]
        {
            sweep_line.debug_order();
            println!("\n{{\"processEvent\": {}}}", event.to_json_debug());
        }
        sorted_events.push(event.clone());

        if operation == Operation::Intersection && event.point.x > rightbound
            || operation == Operation::Difference && event.point.x > sbbox.max.x
        {
            break;
        }

        if event.is_left() {
            let insert_pos = sweep_line.insert(event.clone());
            let insert_pos = match insert_pos {
                Some(insert_pos) => insert_pos,
                None => break
            };

            let index_prev = sweep_line.prev_index(insert_pos);
            let index_next = sweep_line.next_index(insert_pos);
            let maybe_prev = index_prev.map(|idx| sweep_line.get_by_index(idx));
            let maybe_next = index_next.map(|idx| sweep_line.get_by_index(idx));

            compute_fields(&event, maybe_prev, operation);

            let mut mutation_sum = 0;
            if let Some(next) = maybe_next {
                #[cfg(feature = "debug-booleanop")]
                {
                    println!("{{\"seNextEvent\": {}}}", next.to_json_debug());
                }
                let mutations = possible_intersection(&event, &next, event_queue);
                mutation_sum += mutations;
                if mutations == 2 {
                    // Recompute fields for current segment and the one above (in bottom to top order)
                    compute_fields(&event, maybe_prev, operation);
                    compute_fields(&next, Some(&event), operation);
                }
            }

            if let Some(prev) = maybe_prev {
                #[cfg(feature = "debug-booleanop")]
                {
                    println!("{{\"sePrevEvent\": {}}}", prev.to_json_debug());
                }
                let mutations = possible_intersection(&prev, &event, event_queue);
                mutation_sum += mutations;
                if mutations == 2 {
                    let maybe_prev_prev = index_prev
                        .and_then(|idx| sweep_line.prev_index(idx))
                        .map(|idx| sweep_line.get_by_index(idx));
                    // Recompute fields for current segment and the one below (in bottom to top order)
                    compute_fields(&prev, maybe_prev_prev, operation);
                    compute_fields(&event, Some(prev), operation);
                }
            }

            if mutation_sum > 0 {
                sweep_line.fix_rank_range(index_prev.unwrap_or(insert_pos), index_next.unwrap_or(insert_pos));
            }
        } else if let Some(other_event) = event.get_other_event() {
            let index_existing = sweep_line.find(&other_event);
            if let Some(index_existing) = index_existing {
                let index_prev = sweep_line.prev_index(index_existing);
                let index_next = sweep_line.next_index(index_existing);
                let maybe_prev = index_prev.map(|idx| sweep_line.get_by_index(idx));
                let maybe_next = index_next.map(|idx| sweep_line.get_by_index(idx));

                if let (Some(prev), Some(next)) = (maybe_prev, maybe_next) {
                    #[cfg(feature = "debug-booleanop")]
                    {
                        println!("Possible post intersection");
                        println!("{{\"sePostNextEvent\": {}}}", next.to_json_debug());
                        println!("{{\"sePostPrevEvent\": {}}}", prev.to_json_debug());
                    }
                    if possible_intersection(&prev, &next, event_queue) > 0 {
                        sweep_line.fix_rank_range(index_prev.unwrap(), index_next.unwrap());
                    }
                }

                #[cfg(feature = "debug-booleanop")]
                {
                    println!("{{\"removing\": {}}}", other_event.to_json_debug());
                }
                sweep_line.remove_by_index(index_existing);
            } else {
                // This debug assert is only true, if we compare segments in the sweep line
                // based on identity (currently), and not by value (done previously).
                debug_assert!(
                    false,
                    "Sweep line misses event to be removed"
                );
            }
        }
    }

    sorted_events
}
